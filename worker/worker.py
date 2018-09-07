# -*- coding: UTF-8 -*-
import os
import logging

import arrow
import pandas as pd
from bs4 import BeautifulSoup
from toolz.itertoolz import concat
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select

from .data_types import GlycomodComposition, SubmittedMass
from .utils import string_to_dict
from .utils import calc_default_adducts_mono
from .utils import truncated_str_from_dict


class GlycomodWorker:
    """GlycomodWorker accepts string containing N-Glycan masses
    and runs Glycomod search for all specified masses.

    ### So far only MASSES FROM POSITIVE MODE MS CAN BE USED. ###
    ### So far only H+ is USED AS ADDUCT. ###

    Results can be reported as .csv or .txt.
    Masses without matches on Glycomod are reported as NOT FOUND.

    EXAMPLE - from html - all relevant data from web
    'User mass: 1416.64',
    'Adduct ([M+H]+): 1.00727',
    'Derivative mass (Free reducing end): 18.0105546',
    '1397.5080.114(HexNAc)4 (Deoxyhexose)1 (NeuGc)1 (Pent)1',
    '1397.5080.114(Hex)1 (HexNAc)4 (NeuAc)1 (Pent)1',
    '2 structures'

    Data in EXAMPLE gets parsed and stored as objects used for reporting and calculations.

    Monoisotopic or average masses can be used.

    At this point only Da can be used when specifying mass error.
    Default mass error tolerance is 0.5 Da.

    Glycan search properties - monosaccharide presence and number(range) can be
    configured in config.json/"default_mono" and config.json/"default_occurences"
        2 -> "yes" -> monosaccharide must be present
        1 -> "possible" -> monosaccharide may be present
        0 -> "no" -> glycan MUST NOT CONTAIN specified monosaccharide
    BY DEFAULT ONLY N-glycans CONTAINING (Man)3(GlcNAc)2 ARE SELECTED
    """
    def __init__(self, cfg, driver_path="", text_input="", reducing_end=None, adduct="H+", save_txt=False, filename=""):
        self.driver = None
        self.logger = logging.getLogger(name="Worker")
        self.cfg=cfg
        self.driver_path = driver_path
        self.reducing_end_tag = None
        self.reducing_end_mass = 0.0
        self.reducing_end_mass_full = 0.0  # USED AS 'derivative_mass' field input in GM
        self.save_txt = save_txt
        self.filename = filename
        self.text_input = text_input
        self.soup = None
        self.parsed_data = []
        self.compositions = []
        # configure webdriver
        self.driver_options = Options()  
        self.driver_options.add_argument("--headless") 
        try:
            if reducing_end:
                self.reducing_end_tag = reducing_end
                self.reducing_end_mass = self.cfg["reducing_end_tag_mono"][reducing_end]   
                self.reducing_end_mass_full = self.cfg["reducing_end_tag_mono_full"][reducing_end]      
        except ValueError as e:
            self.logger.error(f"End tag {reducing_end} not supported, supporting only 2-AB and ProA at this point\n{e}")
        try:
            self.adduct_info = (adduct, self.cfg["mono_masses_underivatized"][adduct])
        except ValueError as e:
            self.logger.error(f"Adduct {adduct} not supported")

    def output_csv(self):
        out = list(concat([i.prep_csv_out() for i in self.compositions]))
        df = pd.DataFrame.from_records(out, columns=self.cfg["col_names"])
        if self.filename:
            df.to_csv(str(self.filename) + ".csv", index=False)    
        else:
            # set filename in case self.save_txt == True, both files should have the same name
            self.filename = f"results_{arrow.now().format('YYYYMMDD_HH:mm:ss')}"
            df.to_csv(self.filename + ".csv", index=False)
        self.logger.debug(f"Finished saving results as .csv")
    
    def output_text(self, to_std_out=False):
        pretty_txt = self._prettify_text()
        if to_std_out:
            print(pretty_txt)
        else:
            if self.filename:
                with open(str(self.filename) + ".txt", "w") as outfile:
                    outfile.write(pretty_txt)
            else:
                with open(f"results_{arrow.now().format('YYYYMMDD_HH:mm:ss')}.txt", "w") as outfile:
                    outfile.write(pretty_txt)
            self.logger.debug("Finished saving results as .txt")

    def _prettify_text(self):
        """ Transforms raw text to more readable form.
        EXAMPLE RAW:
            'User mass: 1454.0',
            'Adduct ([M+H]+): 1.00727',
            'Derivative mass (Free reducing end): 18.0105546',
            '1434.5020.48(Hex)3 (HexNAc)2 (Deoxyhexose)1 (Pent)3',
            '1435.32-0.337(Hex)1 (HexNAc)3 (Deoxyhexose)2 (Pent)1 (Sulph)3',
            '1435.362-0.379(Hex)4 (HexNAc)1 (Deoxyhexose)2 (Pent)1 (Sulph)2',
            '3 structures'
        EXAMPLE PROCESSED:
            User mass: 1454.0
            Adduct ([M+H]+): 1.00727
            Derivative mass (Free reducing end): 18.0105546
	            1. [MH]+: 1434.502,  Error: 0.48, Comp: (Hex)3(HexNAc)2(Deoxyhexose)1(Pent)3
	            2. [MH]+: 1435.32,  Error: -0.337, Comp: (Hex)1(HexNAc)3(Deoxyhexose)2(Pent)1(Sulph)3
	            3. [MH]+: 1435.362,  Error: -0.379, Comp: (Hex)4(HexNAc)1(Deoxyhexose)2(Pent)1(Sulph)2
            3 structures"""
        prep_text = self.parsed_data.copy()
        for i in prep_text:
            # if unmatched len == 4
            if len(i) > 4:  
                res_string = i[3:len(i)-1]
                prepped = []
                counter = 1
                for res in res_string:                    
                    split_index = res.index("(")
                    num_str, comp_str = res[:split_index:], res[split_index:]
                    comp_str = "".join(comp_str.split())
                    if num_str.find("-") > 0:  # -1 if not found
                        # [mass, error]
                        numbers = [num_str[:num_str.index("-")], num_str[num_str.index("-"):]]
                    else:
                        numbers = [num_str[:num_str.rindex(".")-1], num_str[num_str.rindex(".")-1:]]
                    prepped.append(f"\t{counter}. [MH]+: {numbers[0]:>9},  Error: {numbers[1]:>6}, Comp: {comp_str}")
                    counter += 1
                # replace with processed strings
                i[3:len(i)-1] = prepped
        return "\n\n".join(["\n".join(i) for i in prep_text])

    def _fetch_glycomod_html_data(self):
        """Fetches Glycomod HTML"""
        self.driver = webdriver.Chrome(
            executable_path=self.driver_path, 
            chrome_options=self.driver_options
        )
        self.driver.get(self.cfg["glycomod_link"])
        assert "GlycoMod" in self.driver.title
        self.logger.debug(f"Connected to [{self.driver.title}]")
        # set tolerance
        tol = self.driver.find_element_by_xpath("//input[@name='Tolerance']")
        tol.clear()
        tol.send_keys(self.cfg.get("Tolerance", "0.5"))
        # enter masses into mass textarea
        element = self.driver.find_element_by_xpath("//textarea[@name='Masses']")
        if self.text_input:
            element.send_keys(self.text_input)
        else:
            raise ValueError("FATAL! No data provided for fetching html")
        # select analyte type    
        analyte = Select(self.driver.find_element_by_xpath("//select[@name='Nform']"))
        if self.reducing_end_tag:
            analyte.select_by_visible_text("Derivatised oligosaccharides")
            derivative_name = self.driver.find_element_by_xpath("//input[@name='derivative_name']")
            derivative_name.send_keys(self.reducing_end_tag)
            derivative_mass = self.driver.find_element_by_xpath("//input[@name='derivative_mass']")
            # USING self.reducing_end_mass_full!!!!!!!!!
            derivative_mass.send_keys(str(self.reducing_end_mass_full))         
        else:
            analyte.select_by_visible_text("Free / PNGase released oligosaccharides")
        # select N-glycan monosaccharide content options    
        try:
            for option, value in self.cfg["default_mono"].items():
                if value == 2:
                    Select(
                        self.driver.find_element_by_xpath(f"//select[@name='{option}']")
                        ).select_by_visible_text("yes")
                elif value == 1:
                    Select(
                        self.driver.find_element_by_xpath(f"//select[@name='{option}']")
                        ).select_by_visible_text("possible")
                else:
                    Select(
                        self.driver.find_element_by_xpath(f"//select[@name='{option}']")
                        ).select_by_visible_text("no")
            for option, value in self.cfg["default_occurrences"].items():
                self.driver.find_element_by_xpath(f"//input[@name='{option}']").send_keys(value)
        except Exception as e:
            self.logger.error("Error happened selecting Glycomod options:\n {e}")
            raise

        try:
            self.driver.find_element_by_xpath("//input[@value='Start GlycoMod']").click()
        except Exception as e:
            self.logger.error("Error happended trying to submit Glycomod form:\n{e}")

        self.soup = BeautifulSoup(self.driver.page_source, 'html5lib')
        self.driver.close()
    
    def _parse_gm_html(self) -> list: 
        """Parses HTML for relevant data about glycan compositions"""
        if self.soup:
            items = []
            glycans_list = []
            for hr in self.soup.find_all("hr"):
                for item in hr.find_next_siblings():
                    if item.name == 'hr':
                        break
                    items.append(item.text)
            # this is a bit silly but I do not want to use regexes
            clean_text = ''.join(items).replace(u'\nglycoform mass\nΔmass (Dalton)\nstructure\ntype\nLinks', "")\
                .replace("high_manUniCarbKB", "").replace("hybrid/complexUniCarbKB", "").\
                replace(" -UniCarbKB", "").replace("high_man", "").replace("hybrid/complex", "").\
                replace("paucimannose", "").\
                replace(" -", "").\
                replace("\n\n\n\nSIB Swiss Institute of Bioinformatics | Disclaimer", "").\
                replace("Back to the Top\n\n", "")
            unprocessed_strings = clean_text.split(' found.', clean_text.count('found.')-1)
            for i in unprocessed_strings:
                ll = i.split("\n")
                ll2 = []
                for j in range(len(ll)):
                    if len(ll[j]) > 1:
                        ll2.append(ll[j])
                glycans_list.append(ll2)           
            self.parsed_data = glycans_list 
        else:
            raise ValueError("NO GLYCOMOD DATA WAS FETCHED.")
    
    def _create_glycan_objects(self):
        """Creates SubmittedMass objects from parsed html data"""
        for result in self.parsed_data:
            submitted = 0.0
            compositions = []
            for element in result:
                if element.startswith('User mass: '):
                    submitted = element[len('User mass: '):]
                    if " " in submitted:
                        submitted = "".join(submitted.split())

                if element[0].isdigit():
                    if "0 structures" in element:
                        compositions.append(
                            GlycomodComposition(
                                theoretical_MH=0.0,
                                delta=1000.0,
                                long_notation="NOT FOUND",
                                short_notation="NOT FOUND",
                                theoretical_MTagH=0.0,
                                theoretical_MTagNa=0.0,
                                theoretical_MTagK=0.0,
                                theoretical_MTagH2=0.0,
                                theoretical_MTagHNa=0.0,
                                theoretical_MTagHK=0.0,
                                theoretical_MTagNa2=0.0,
                                theoretical_MTagNH4=0.0,
                            )
                        )
                    elif "structure" not in element:
                        split_index = element.index("(")
                        num_str, comp_str = element[:split_index:], element[split_index:]
                        # find returns -1 if not found
                        if num_str.find("-") > 0:  
                            numbers = [num_str[:num_str.index("-")], num_str[num_str.index("-"):]]
                        else:
                            # if positive, python handles .17 as 0.17, so no problem
                            numbers = [num_str[:num_str.rindex(".")], num_str[num_str.rindex("."):]]
                        comp_dict = string_to_dict(comp_str)
                        # if self.use_avg_vals: # TODO
                        #    raise NotImplementedError
                        adduct_ions = calc_default_adducts_mono(comp_dict, self.cfg, reducing_end=self.reducing_end_tag)
                        compositions.append(
                            GlycomodComposition(
                                theoretical_MH=float(numbers[0]),
                                theoretical_MTagH=adduct_ions["H+"],
                                theoretical_MTagNa=adduct_ions["Na+"],
                                theoretical_MTagK=adduct_ions["K+"],
                                theoretical_MTagH2=adduct_ions["2H2+"],
                                theoretical_MTagHNa=adduct_ions["HNa2+"],
                                theoretical_MTagHK=adduct_ions["HK2+"],
                                theoretical_MTagNa2=adduct_ions["2Na2+"],
                                theoretical_MTagNH4=adduct_ions["NH4+"],
                                delta=float(numbers[1]),
                                long_notation=comp_str,
                                short_notation=truncated_str_from_dict(comp_dict)
                            )
                        )
            self.compositions.append(
                SubmittedMass(
                    experimental_mass=submitted,
                    # peak_number=peak_num,
                    adduct=self.adduct_info[0],
                    adduct_mass=self.adduct_info[1],
                    red_end_tag=self.reducing_end_tag,
                    red_end_tag_mass=self.reducing_end_mass,
                    glycomod_structures=compositions
                )
            )

    def run(self):
        """Run GlycomodWorker search and report results"""
        self._fetch_glycomod_html_data()
        self._parse_gm_html()
        self._create_glycan_objects()
        self.output_csv()
        if self.save_txt:
            try:
                self.output_text()
            except Exception as e:
                self.logger.error("Error outputing text: {e}")

if __name__ == "__main__":
    from pprint import pprint
    import json
    with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.json")), "r") as cc:
        conf = json.load(cc)
    gw = GlycomodWorker(text_input="911.30\n1057.33", cfg=conf, reducing_end=None) 
    gw.run()
    gw.output_text(to_std_out=True)
    # gw._fetch_glycomod_html_data()
    # with open("./tests/test_html/minimal_test_w_multi_matched.html", "r") as ff:
    #     kff = ff.read()
    # soup = BeautifulSoup(kff, 'html5lib')
    # gw.soup = soup
    # gw._parse_gm_html()
    # gw._create_glycan_objects()
    # # pprint(gw.compositions)
    # gw.output_csv()
    # gw.output_text(to_std_out=True)