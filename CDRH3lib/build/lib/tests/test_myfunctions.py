from qualiloop import myfunctions
import pandas as pd
import os


def test_get_input_seq():
    current=os.getcwd()
    file=open(os.path.join(current,"test_PDB_1F11_1.pdb.model"))
    to_test=myfunctions.get_input_seq(file)
    correct_test=pd.read_csv(os.path.join(current,"test_get_input_seq_correct.csv"), header=0)
    assert to_test == correct_test


def test_get_loop_seq():
	correct_test=[["95","D"],["96","Y"],["97","G"],["98","S"],["99","T"],["100","Y"],
	["100A","G"],["100B","F"],["101","A"],["102","Y"]]
	

