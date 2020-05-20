# Kernel Methods for High-Resolution Genomic Prediction across Multiple Environments
Core of functions to build gaussian kernel, arc-cosine and GBLUP with additive, dominance effects and considering environmental information

# Contents
   * [1. Overview](#p0)
   * [2. Data sets](#p2)
      * [2.1 Helix Seeds](#p2.1)
      * [2.2 University of São Paulo (USP)](#p2.2)
   * [3.  Environmental Typing](#p3)
      * [3.1  Collecting data](#p3.1)
      * [3.2  Matrix of Environmental Covariables](#p3.1)
   * [4. Kernel Methods](#p4)
       * [4.1  GBLUP](#p4.1)
       * [4.2  Gaussian Kernel (GK)](#p4.2)
       * [4.3  Deep Kernel Kernel (DK)](#p4.3)      
   * [5. Statistical Models](#p5)
   * [6. Genomic Prediction](#p6)
   * [6. Variance Components](#p7)
   * [7. References](#p8)
   * [8. Authorship and Acknowledgments](#p9)
 
 
 ----------------------------------------------------------------------------------------------------------------
<div id="p0" />

# Overview

This page was developed to provide genomic prediction codes using environmental data, dominance effects and reaction norms. We also provide codes for building kernels such as GBLUP, Gaussian Kernel and Deep Kernel. The base article for these codes is the paper of Costa-Neto et al. (submitted to Heredity journal in May 2020). All codes are free and any questions, suggestions or criticism to improve the codes please look for the authors of this page. Thank you very much!

 ----------------------------------------------------------------------------------------------------------------
<div id="p2" />

# Data sets

<div id="p2.1" />

## Helix Seeds

> 452 tropical maize single-crosses provided by Helix Sementes®, Brazil. Hybrids were obtained from crosses between 128 inbred lines and were evaluated for grain yield (GY) and plant height (PH). Field trials were carried out using a randomized complete block design with two replicates each, allocated across five sites for GY and three for PH during the growing season of 2015.

> Inbred lines were genotyped via the Affymetrix® Axiom® Maize Genotyping Array (Unterseer et al. 2014) with 660K SNP markers. Quality control for SNPs was made based on call rate (CR), in which all markers with any missing data point were excluded, and minor allele frequency (MAF) procedures, in which markers with a low level of polymorphism (MAF < 0.05) were removed. Hybrid genotypes were scored by an allelic combination of homozygous markers of parental lines. After quality control, 37,625 SNP were used to compare the imputation methods.

Check the [CIMMYT DATA VERSE Repository](https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10887) for the full genotypic and phenotypic data


<div id="p2.2" />

## USP data set

> 906 maize single-crosses obtained from a full dial- lel, according to Griffing’s method 4, divided into two heterotic groups, flint and dent, with 34 and 15 lines, respec- tively. Moreover, each heterotic group has a representative line, frequently used as the tester in our breeding program.

> The experimental scheme used to evaluate the hybrids was an augmented block design (unreplicated trial) consisted of small blocks, each with 16 unique hybrids and two checks. Trials were carried out in Anhembi (22°50′51′′S, 48°01′06′′W, 466 m) and Piracicaba, at São Paulo State, Brazil (22°42′23′′S, 47°38′14′′W, 535 m), during the second growing season of 2016 and 2017, cultivated between January to June. In both sites and years, the hybrids were evaluated under two nitrogen (N) levels, low (LN) with 30 kg N ha−1, and normal (NN) with 100 kg N ha−1.

> The genotyping of the 49 tropical inbred lines was per- formed by Affymetrix® platform, containing about 614,000 SNPs (Unterseer et al. 2014). Then, markers with low call rate (< 95%), minor allele frequency (MAF < 0.05) and heterozygous loci on at least one individual were removed. The missing markers were imputed using the [snpReady](https://github.com/italo-granato/snpReady) R package. Finally, the resulting 146,365 SNPs high-quality polymorphic SNPs were used to build the artificial hybrids genomic matrix, deduced by combining the genotypes from its two parents.

Check the [Mendeley Repository](https://data.mendeley.com/datasets/tpcw383fkm/3) for the full genotypic and phenotypic data


--------------------------------------------------------------------------------------------------------------------------------
<div id="p3" />

# Environmental Typing

The environmental typing (envirotyping) pipeline were conducted using the functions of the [EnvRtype](https://github.com/allogamous/EnvRtype) R package. This package has functions for supporting the collection of environmental data **get_weather()**, processing environmental data **processWTH()** and build of the **W** matrix of envirotype covariables **W.matrix()**. Finally, this package helps the construction of the genomic x envirotyping kinships using **get_kernels()**  functions, which easily can run into a Bayesian Genomic-enabled Prediction framework implemented in [BGGE](https://github.com/italo-granato/BGGE) package. Bellow we present a brief example of the use of both packages to run a genomic prediction considering reaction norms.

<div id="p3.1" />

## Collecting data

> Environmental data were obtained from NASA POWER data base using the function  **get_weather()** based on the geographic coordinates and planting dates for each environment:


<div id="p3.2" />

## Matrix of Environmental Covariables

> Additional environmental variables describing ecophysiological processes (e.g., evapotranspiration, effect of temperature on radiation use efficiency) were computed using the function **processWTH()**:

> Then, it's possible to create an environmental covariable matrix **W**, with q environments and k combinations of time intervals (e.g., phenologycal stages) x quantiles (distribution of the environmental data) x environmental factor, as follows:


 ----------------------------------------------------------------------------------------------------------------
<div id="p4" />

# Kernel Methods

Below are the codes for three types of kernel methods. We use these methods to model additive (**A**), dominance (**D**) and environmental (**W**) effects. On this page we will make the codes available, especially for obtaining the **D** effects, but we will exemplify the kernels using only the **A** effects. To run **D** and **W** just replace the matrix A with the respective **D** and **W**.

> Obtaining dominance effects is given by:


<div id="p4.1" />

## GBLUP

> Relationship Kernels based on the Genomic Best Linear Unbiased Predictior (GBLUP) can be implemented as:


<div id="p4.2" />

## Gaussian Kernel (GK)

> Relationship Kernels based on Gaussian Kernel (GK) can be implemented as:

> If the user want to compute the bandwith factor, a likehood marginal function can be implemented by providing the phenotypic data of the training set as:


<div id="p4.3" />

## Deep Kernel (DK)

> Relationship Kernels based on Deep Kernels (DK) cab be implemented based on the arc-cosine method presented in genomic prediction by Cuevas et al (2019) and Crossa et al (2020). Firstly, a base arc-cosine kernel are computed using the molecular matrix data (coded as additive = (0,1,2) or dominance) or environmental data (per environment or per genotype-environment combinations) as follows:

> Then, using the function **get_** we can


 ----------------------------------------------------------------------------------------------------------------
<div id="p5" />

# Statistical Models

Five genomic prediction models were presented using the function **get_kernels** from EnvRtype package.

### Main Additive Effect Model (without GE effects)


### Main Additive plus Dominance Effects Model (without GE effects)

### Main Additive-Dominance effects plus GE deviation (GE = AE + DE)


### Main Additive-Dominance effects plus Envirotyping information (W)

### Main Additive-Dominance effects plus GE reaction norm (W+AW+DW)


 ----------------------------------------------------------------------------------------------------------------
<div id="p6" />

## Genomic Prediction

> Genomic predictions were performed using the Bayesian Genotype plus Genotype × Environment (BGGE) package (Granato et al. 2018) fitted to 10,000 iterations with the first 1,000 cycles removed as burn-in with thinning equal to 2. 


 ----------------------------------------------------------------------------------------------------------------
<div id="p7" />

## Variance Components

> Variance components can be extracted using the function **Vcomp.BGGE** provided in this [link]() or using the following source code:



 ----------------------------------------------------------------------------------------------------------------
<div id="p8" />

# References

* Crossa J, Martini JWR, Gianola D, et al (2019) Deep Kernel and Deep Learning for Genome-Based Prediction of Single Traits in Multienvironment Breeding Trials. Front Genet 10:1–13. https://doi.org/10.3389/fgene.2019.01168

* Cuevas J, Montesinos-López O, Juliana P, et al (2019) Deep Kernel for genomic and near infrared predictions in multi-environment breeding trials. G3 Genes, Genomes, Genet 9:2913–2924. https://doi.org/10.1534/g3.119.400493

* Granato I, Cuevas J, Luna-vázquez F, et al (2018) BGGE : A New Package for Genomic-Enabled Prediction Incorporating Genotype x Environment Interaction Models. 8:3039–3047. https://doi.org/10.1534/g3.118.200435

* Pérez-Elizalde S, Cuevas J, Pérez-Rodríguez P, Crossa J (2015) Selection of the Bandwidth Parameter in a Bayesian Kernel Regression Model for Genomic-Enabled Prediction. J Agric Biol Environ Stat 20:512–532. https://doi.org/10.1007/s13253-015-0229-y

* Pérez-Rodríguez P, Gianola D, González-Camacho JM, et al (2012) Comparison between linear and non-parametric regression models for genome-enabled prediction in wheat. G3 Genes, Genomes, Genet. https://doi.org/10.1534/g3.112.003665




 ----------------------------------------------------------------------------------------------------------------
<div id="p9" />

# Authorship and Acknowledgments

For any questions, suggestions or corrections, please contact the author of this tutorial at **Germano Costa neto** germano.cneto@usp.br

We thank prof. Jaime Cuevas for his availability for Deep Kernel codes. For more details on this subject, contact him or look for more information in Cuevas et al. (2019).

We thank CIMMYT and the University of São Paulo (USP) for their support in the development of the project.

The theories exposed on this page are presented in the work of Costa Neto et al (not published yet) and were developed under the authorship of Germano Costa Neto, Jose Crossa and Roberto Fritsche Neto.
    


