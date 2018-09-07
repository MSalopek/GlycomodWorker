# -*- coding: UTF-8 -*-
import re
from http.client import HTTPConnection
from typing import Dict
import hashlib

# Cytonize all util functions??

def string_to_dict(string)-> Dict[str, int]:  
    elements_count = {}
    # this should never happen
    if string == "(Man)3(GlcNAc)2":  
        elements_count["(Man)"] = 3
        elements_count["(GlcNAc)"] = 2
        return elements_count

    elements = string.split()
    if "+" in string:
        plus_split = string.split(" + ")
        split_elems = plus_split[0].split()
        for i in split_elems:
            sym = i.strip("0123456789")
            num = i[len(sym):]
            elements_count[sym] = int(num)
        plus_split[1] = plus_split[1].rstrip()
        if plus_split[1] == "(Man)3(GlcNAc)2":  
            # USUALLY there's only "(Man)3(GlcNAc)2" after " + "
            elements_count["(Man)"] = 3
            elements_count["(GlcNAc)"] = 2
        else:
            raise ValueError(f"Unexpected trailing sequence {plus_split[1]}, expected (Man)3(GlcNAc)2")
        return elements_count 

    for i in elements:
        sym = i.strip("0123456789")
        num = i[len(sym):]
        elements_count[sym] = int(num)
    return elements_count


def truncate_str(string) -> str:
    """Takes string representing long form glycan composition (from Glycomod) and returns short form composition
        ### EXAMPLE ###
        >>>truncate_str("(Hex)10 (HexNAc)5 (NeuAc)5 + (Man)3(GlcNAc)2")) 
        >>>"H13N7S5"
    """
    str_dict = string_to_dict(string)
    total_hex = str_dict.get("(Man)", 0) + str_dict.get("(Hex)", 0)
    total_hexnac = str_dict.get("(HexNAc)", 0) + str_dict.get("(GlcNAc)", 0)
    return f"{'H' + str(total_hex)}{'N' +str(total_hexnac)}" + \
        f"{'F' + str(str_dict['(Deoxyhexose)']) if str_dict.get('(Deoxyhexose)') else ''}" + \
        f"{'P' + str(str_dict['(Pent)']) if str_dict.get('(Pent)') else ''}" + \
        f"{'S' +  str(str_dict['(NeuAc)']) if str_dict.get('(NeuAc)') else ''}" + \
        f"{'Sg' +  str(str_dict['(NeuGc)']) if str_dict.get('(NeuGc)') else ''}" + \
        f"{'KDN' + str(str_dict['(KDN)']) if str_dict.get('(KDN)') else ''}"+ \
        f"{'Sulph' + str(str_dict['(Sulph)']) if str_dict.get('(Sulph)') else ''}" + \
        f"{'Phos' + str(str_dict['(Phos)']) if str_dict.get('(Phos)') else ''}" + \
        f"{'HexA' + str(str_dict['(HexA)']) if str_dict.get('(HexA)') else ''}"


def truncated_str_from_dict(dictionary):
    total_hex = dictionary.get("(Man)", 0) + dictionary.get("(Hex)", 0)
    total_hexnac = dictionary.get("(HexNAc)", 0) + dictionary.get("(GlcNAc)", 0)
    return f"{'H' + str(total_hex)}{'N' +str(total_hexnac)}" + \
        f"{'F' + str(dictionary['(Deoxyhexose)']) if dictionary.get('(Deoxyhexose)') else ''}" + \
        f"{'P' + str(dictionary['(Pent)']) if dictionary.get('(Pent)') else ''}" + \
        f"{'S' +  str(dictionary['(NeuAc)']) if dictionary.get('(NeuAc)') else ''}" + \
        f"{'Sg' +  str(dictionary['(NeuGc)']) if dictionary.get('(NeuGc)') else ''}" + \
        f"{'KDN' + str(dictionary['(KDN)']) if dictionary.get('(KDN)') else ''}"+ \
        f"{'Sulph' + str(dictionary['(Sulph)']) if dictionary.get('(Sulph)') else ''}" + \
        f"{'Phos' + str(dictionary['(Phos)']) if dictionary.get('(Phos)') else ''}" + \
        f"{'HexA' + str(dictionary['(HexA)']) if dictionary.get('(HexA)') else ''}" 


def calc_theor_mono_mass(dictionary, cfg, prec=6, reducing_end=None) -> float:
    """Returns theoretical monoisotopic mass for glycan in dictionary form
        - reducing end type specified by reducing_end parameter ("2-AB" or "ProA") ELSE treated as unlabeled nonreduced"""
    reducing_end_tag_mass = 0.0
    if reducing_end is not None:
        if reducing_end in cfg["reducing_end_tag_mono"].keys():
            reducing_end_tag_mass = cfg["reducing_end_tag_mono"][reducing_end]
        else:
            raise ValueError(f"Invalid reducing end tag. Available tags:  {list(cfg['reducing_end_tag'].keys())}")
    # I don't know why there is H3O+, but it is in accordance with what Glycomod and Glycoworkbench show as mono masses
    return round(
        sum(
            [
                dictionary[key]*cfg["mono_masses_underivatized"][key]
                for key in dictionary.keys()
            ]
        ) + cfg["mono_masses_underivatized"]["H3O+"] + reducing_end_tag_mass,
        prec
    )


def _calc_theor_mono_mass_adducts(dictionary, cfg, charge, reducing_end_mass, adduct_mass, prec=6):
    return round(
        (sum(
            [
                dictionary[key]*cfg["mono_masses_underivatized"][key]
                for key in dictionary.keys()
            ]
        ) + cfg["mono_masses_underivatized"]["H2O"] + reducing_end_mass + adduct_mass) / charge,
        prec
    )


def calc_default_adducts_mono(dictionary, cfg, prec=4, reducing_end=None):
    reducing_end_tag_mass = 0.0
    if reducing_end is not None:
        if reducing_end in cfg["reducing_end_tag_mono"].keys():
            reducing_end_tag_mass = cfg["reducing_end_tag_mono"][reducing_end]
        else:
            raise ValueError(f"Invalid reducing end tag. Available tags: {list(cfg['reducing_end_tag'].keys())}")
    adduct_ions = {
        i: _calc_theor_mono_mass_adducts(dictionary, cfg, j[1], reducing_end_tag_mass, j[0], prec=prec) for i,j in cfg["adducts"].items()
    }
    return adduct_ions


def calc_theor_avg_mass(dictionary, cfg, prec=6, reducing_end=None) -> float:
    """Returns theoretical average mass for glycan in dictionary form"""
    reducing_end_tag_mass = 0.0
    if reducing_end is not None:
        if reducing_end in cfg["reducing_end_tag_avg"].keys():
            reducing_end_tag_mass = cfg["reducing_end_tag_avg"][reducing_end]
    return round(
        sum(
            [
                dictionary[key]*cfg["avg_masses_underivatized"][key] 
                for key in dictionary.keys()
                ]
            ) + cfg["avg_masses_underivatized"]["H3O+"] + reducing_end_tag_mass,
        prec
    )


def check_internet_conn():
    conn = HTTPConnection("www.google.com", timeout=5)
    try:
        conn.request("HEAD", "/")
        conn.close()
        return True
    except:
        conn.close()
        return False


def validate_filename(string):
    filename = ""
    if len(string) > 1:
        for ind, char in enumerate(string):
            if char.isalnum():
                filename = filename + char
            elif char == " " or char == "-" or char == "_":
                filename = filename + char
            else:
                raise ValueError(f"Char {ind} invalid in {string}. File name must consist of letters, numbers '-' and '_'. Please provide a valid file name")
    return filename


def check_text_input_consistency(string):
    stripped_text = re.sub(r'\s+', '', string)
    for char in stripped_text:
        if char.isdigit() or char == ".":
            continue
        else:
            raise ValueError(f"Your input contains invalid character '{char}'; Input should consist of floating point numbers.")
