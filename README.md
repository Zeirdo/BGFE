# BGFE

A Deep Learning Model for ncRNA-Protein Interaction Predictions Based on Improved Sequence Information

we proposed a computational method using the deep learning stacked autoencoder network algorithm to mine the hidden high-level features from ncRNA and protein sequences and fed them into random forest predictors to predict ncRNA binding proteins, stacked assembling is used to improve performance. The 5-fold cross-validation and widely used evaluation measure are used to evaluate our method's performance.

## Dependency: 

blastpgp: to convert a protein sequence into a PSSM matrix;

Matlab: to extract Zernike moment features;

deep learning lib keras: https://github.com/fchollet/keras/  (version Keras-0.1.2)

    #If something goes wrong, this is because keras's API has changed a lot since the 1.0 upgrade

    #you need can able to modify the contents of .../keras/layers/convolutional.py

    #or you can ask for it by e-mail

machine learning lib scikit-learn: https://github.com/scikit-learn/scikit-learn 

MinGW and libpython are also required,if in the Windows environment

## Usage:

python BGFE.py -dataset=RPI2241  

BGFE will do 5-fold cross-validation for it. you can also choose other datasets, such as RPI1807, RPI488, RPI2241. 

python BGFE.py -r=RNA_fasta_file -p=protein_fasta_file   

it will predict pairwise interaction score for RNAs and protiens in input file.


Contact: TS16170022A3@cumt.edu.cn
