# GPMsDB-dbtk

GPMsDB-dbtk v1.0.1 was released on June xx, 2022. 

GPMsDB-tk/GPMsDB-dbtk are software toolkits for assigning taxonomic identification to user-provided MALDI-TOF mass spectrometry profiles obtained from bacterial and archaeal cultured isolates. They take advantages of a newly developed database of protein mass profiles predicted from ~200,000 bacterial and archaeal genome sequences. This toolkit is also designed to work with customized databases, allowing microbial identification based on user-provided genome/metagenome-assembled genome (MAG) sequences. The ms-identification-tk is open source and released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. 

Please post questions and issues related to ms-identification-tk on the Issues section of the GitHub repository.

## Installing and using GPMsDB-dbtk

Prerequisites
* Python (version 3.7 or higher)
* Cython (version 0.29.1 or higher)
* biolib (https://github.com/dparks1134/biolib) (version 0.1.8 or higher)
* matplotlib (version 3.5.0 or higher)
* Prodigal (https://github.com/hyattpd/Prodigal) (version 2.6.3 or higher)
* HMMER (https://github.com/EddyRivasLab/hmmer) (version 3.3.2 or higher)

In the source directory, the following command will compile and install the software in your python environment.
```bash
python setup.py install
```

GPMsDB-dbtk requires ~60 GB of external data that needs to be downloaded and unarchived:
```bash
wget https://....data.tar.gz
tar xvzf ....data.tar.gz
```

GPMsDB-dbtk requires an environment variable named GPMsDB_PATH to be set to the directory containing the unarchived reference data.
```bash
export GPMsDB_PATH=/path/to/release/package/
```

### Features

* Genome(s) to massDB:
  * genome_wf     -> Full genomes to ms data workflow
  * list_db       -> List genome entries in the custom ms database
  * update_db     -> Add peak_list(s) to the custom ms database
  * remove_genome -> Delete entries from the custom ms database
			
## Bug Reports

Please report bugs through the GitHub issues system, or contact Yuji Sekiguchi (y.sekiguchi@aist.go.jp)

## Copyright

Copyright (C) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)

This package (the majority of the scripts in the package) is under the conditions of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. See LICENSE for further details.