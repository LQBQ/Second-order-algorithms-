# EEM Fluorescence and Chemometrics for Olive Oil Adulteration Detection

# About the Project

This repository contains supplementary material for the work "EEM fluorescence coupled with second-order algorithms in R for detecting extra virgin olive oil adulteration: an experimental approach for undergraduate students - Part II", submitted to the Journal of Chemical Education. The study aims to present an experimental approach to teaching second-order chemometric algorithms in the context of detecting adulteration in extra virgin olive oil (EVOO).

# Motivation

The adulteration of extra virgin olive oil is a recurring problem in the food industry, requiring advanced analytical methods for efficient detection. This study proposes the use of excitation-emission fluorescence spectroscopy (EEM) combined with chemometric algorithms such as PARAFAC, Tucker-3, and nPLS to identify and quantify adulterants in EVOO. The main differentiator of this work is its educational approach, aimed at undergraduate chemistry students, enabling the understanding and practical application of these methodologies.

# Institution and Team

This project was developed at the Laboratory of Biological Chemistry and Chemometrics (LQBQ), part of the Institute of Chemistry at the Federal University of Rio Grande do Norte (UFRN), under the coordination of Prof. Dr. Kássio M.G. de Lima.

# Authors:

Ramon B. P. Alves (LQBQ - UFRN)

Anne B. F. Câmara (LQBQ - UFRN)

Hellyda ... (LQBQ - UFRN)

Ayrton ... (LQBQ - UFRN / Laboratory...)

Luiz S. Neves (Institute of Chemistry - UFRN)

Kássio M.G. de Lima (LQBQ - UFRN)

Contact: kassio.lima@ufrn.br

# Methodology

Samples: 15 adulterated EVOO samples were prepared with soybean oil at concentrations ranging from 0.5% to 10% (w/w) and 12 pure EVOO samples.

Data Acquisition: EEM matrices were obtained using a Shimadzu RF-5301PC spectrofluorometer.

Chemometric Analysis: PARAFAC, Tucker-3, and nPLS models were built in R (version 4.3.1) using packages such as multiaway, ThreeWay, pls, among others.

Validation: The models were evaluated using metrics such as the coefficient of determination (R²), root mean square error (RMSEC/RMSEP), and classification accuracy.

# Available Resources in the Repository

This repository contains:

R Scripts for importing, processing, and modeling data.

Tutorials for installing R and the necessary packages.

Database containing fluorescence matrices of the samples.

Reports and supplementary materials on the applied models and their results.

Questionnaire applied to students for evaluating the educational approach.

# Educational Application

This approach was implemented in a Chemometrics course for undergraduate Chemistry students at UFRN. The students were introduced to the use of EEM and second-order chemometric algorithms, replicating the models with a real colorectal cancer dataset. The evaluation showed positive feedback, indicating that the methodology contributed to the students' academic and professional training.

# How to Use

Install R and the required packages listed in the tutorial.

Import the data from the fluorescence matrix database.

Run the scripts to process the data and build the chemometric models.

Analyze the results using the provided metrics.

# Final Considerations

This project aims not only to contribute to food quality control but also to expand the teaching of chemometrics at the higher education level. The use of R as a free and accessible tool strengthens this initiative, making chemometric methods more widespread among students and researchers.

# References

The full list of references used in this study can be found in the associated article and the supplementary materials of this repository.

We hope this material will be useful for researchers and educators interested in applying chemometrics to real analytical problems!

For questions or suggestions, please contact us via email.

License: This repository is available under the MIT license. Feel free to contribute, modify, and share.
