# What is GlycomodWorker?

GlycomodWorker is a tool for fetching data from Glycomod(ExPASy). It avoids all the copy pasting when using Glycomod for composition search with experimental mass spectrometry data. Results are stored as .csv tables with relevant adduct info calculated for you.

## Input 
Text files (.txt) containing newline separated experimental masses. All masses must be floating point numbers.

## CLI usage
Navigate to a directory where GlycomodWorker is saved.
```sh
curdir:$ python -m GlycomodWorker --help
```
To check available commands.

-------------
BASIC USE:
```sh
curdir:$ python -m GlycomodWorker /home/ms/GlycomodWorkflow/GlycomodWorker/test.txt --tag ProA --echo
```
By default results are stored in the directory where GlycomodWorker is located.
File name example:
`results_20180512:07:31:55.csv`