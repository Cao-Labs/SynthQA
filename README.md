# SynthQA - Hierarchical Machine Learning-basedProtein Quality Assessment
## Dependancies
This tool requires Python 3.7, as well as the latest versions of the following:
- [ANDIS (standalone version)](http://qbp.hzau.edu.cn/ANDIS/)  
- [PyRosetta4](http://www.pyrosetta.org/dow)
- Numpy
- Scipy
- Scikit-learn
- h5py

## Installation
**NOTE**: This tool was tested on `Ubuntu 20.04 LTS`. Installation and requirements may vary on other operating systems.

After installing ANDIS, be sure to modify execution permissions with this
command: `chmod +x ANDIS`.

In main.py, please change the following variables found at the top of the file:
- Change **QA_DIR** to the absolute path of this tool's directory.
- Change **ANDIS_PATH** to the absolute path of your installed ANDIS directory.

You can install the required Python libraries using requirements.txt:
`pip3 install -r requirements.txt`

Remember to also download and install PyRosetta4 (http://www.pyrosetta.org/dow) before proceeding.

## Running
Example usage:
```
python main.py -i ./input_folder_name -o output_folder_name
```
The input folder must contain at least one protein model (PDB format) for this tool to work.

## Testing
A sample output can be found in `example/sample_output`. It was generated using the following command:
```
python main.py -i ./example/T0953s2 -o ./example/sample_output
```


--------------------------------------------------------------------------------------
Developed by Mikhail Korovnik and Prof. Renzhi Cao at Pacific Lutheran University:

Please contact Renzhi Cao for any questions: caora@plu.edu (PI)
