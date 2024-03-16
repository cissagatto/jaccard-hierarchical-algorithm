# Jaccard

This code generate partitions for a multilabel dataset using the Jaccard Index similarity measure. We use HCLUST with 6 linkage metrics to generate several partitions. You may build the partition with the highest coefficient. This code also provide an analysis about the partitioning

## Preparing your experiment

### Step-1
Confirms if the folder .......

### Step-2
Copy this code and place it where you want. The folder configurations is "~/jaccard"

### Step-3
A file called _datasets-original.csv_ must be in the *root project folder*. This file is used to read information about the datasets and they are used in the code. We have 90 multilabel datasets in this _.csv_ file. If you want to use another dataset, please, add the following information about the dataset in the file:


| Parameter    | Status    | Description                                           |
|------------- |-----------|-------------------------------------------------------|
| Id           | mandatory | Integer number to identify the dataset                |
| Name         | mandatory | Dataset name (please follow the benchmark)            |
| Domain       | optional  | Dataset domain                                        |
| Instances    | mandatory | Total number of dataset instances                     |
| Attributes   | mandatory | Total number of dataset attributes                    |
| Labels       | mandatory | Total number of labels in the label space             |
| Inputs       | mandatory | Total number of dataset input attributes              |
| Cardinality  | optional  |                                                       |
| Density      | optional  |                                                       |
| Labelsets    | optional  |                                                       |
| Single       | optional  |                                                       |
| Max.freq     | optional  |                                                       |
| Mean.IR      | optional  |                                                       | 
| Scumble      | optional  |                                                       | 
| TCS          | optional  |                                                       | 
| AttStart     | mandatory | Column number where the attribute space begins*       | 
| AttEnd       | mandatory | Column number where the attribute space ends          |
| LabelStart   | mandatory | Column number where the label space begins            |
| LabelEnd     | mandatory | Column number where the label space ends              |
| Distinct     | optional  |                                                       |
| xn           | mandatory | Value for Dimension X of the Kohonen map              | 
| yn           | mandatory | Value for Dimension Y of the Kohonen map              |
| gridn        | mandatory | X times Y value. Kohonen's map must be square         |
| max.neigbors | mandatory | The maximum number of neighbors is given by LABELS -1 |

* Because it is the first column the number is always 1.


## STEP 4
You need to have installed all the R packages required to execute this code on your machine. Check out which are needed in the file *libraries.R*. This code does not provide any type of automatic package installation! You can use the Conda environment that I created to perform this experiment. Below are the links to download the files.

| [download txt](https://www.4shared.com/s/fUCVTl13zea) | [download yml](https://www.4shared.com/s/f8nOZyxj9iq) | [download yaml](https://www.4shared.com/s/fk5Io4faLiq) |


# Run

To run, first enter the folder ~/jaccard/R in a terminal and the type: (check this information)

```
Rscript jaccard.R [number_dataset] [number_cores] [number_folds] [validation] [folder]
```

Where:

_number_dataset_ is the dataset number in the datasets.csv file

_number_cores_ is the number of cores that you wanto to use in paralel

_number_folds_ is the number of folds you want for cross-validation

_validation_ 0 if you dont want the validation set and 1 if you want

_folder_ temporary folder like SHM or SCRATCH to speed up the process

## Folder Structure



## Acknowledgment
- This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001.
- This study was financed in part by the Conselho Nacional de Desenvolvimento Científico e Tecnológico - Brasil (CNPQ) - Process number 200371/2022-3.
- The authors also thank the Brazilian research agencies FAPESP financial support.

# Contact
elainececiliagatto@gmail.com

## Links

| [Site](https://sites.google.com/view/professor-cissa-gatto) | [Post-Graduate Program in Computer Science](http://ppgcc.dc.ufscar.br/pt-br) | [Computer Department](https://site.dc.ufscar.br/) |  [Biomal](http://www.biomal.ufscar.br/) | [CNPQ](https://www.gov.br/cnpq/pt-br) | [Ku Leuven](https://kulak.kuleuven.be/) | [Embarcados](https://www.embarcados.com.br/author/cissa/) | [Read Prensa](https://prensa.li/@cissa.gatto/) | [Linkedin Company](https://www.linkedin.com/company/27241216) | [Linkedin Profile](https://www.linkedin.com/in/elainececiliagatto/) | [Instagram](https://www.instagram.com/cissagatto) | [Facebook](https://www.facebook.com/cissagatto) | [Twitter](https://twitter.com/cissagatto) | [Twitch](https://www.twitch.tv/cissagatto) | [Youtube](https://www.youtube.com/CissaGatto) |

# Thanks
