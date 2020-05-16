# Kernel Methods for High-Resolution Genomic Prediction across Multiple Environments
Core of functions to build gaussian kernel, arc-cosine and GBLUP with additive, dominance effects and considering environmental information

# Contents

   * [1. HEL data set](#p1)
   * [2. USP data set](#p2)
   * [3. Environmental Typing](#p3)
   * [4. Kernel Methods](#p4)
   * [5. Statistical Models](#p5)
   * [6. Cross-Validation Schemes](#p6)
   * [7. A Toy Example](#p7)
   * [8. References](#P8)
   * [9. Authorship and Acknowledgments](#P8)
 
 ----------------------------------------------------------------------------------------------------------------
<div id="p1" />

# HEL data set

452 tropical maize single-crosses provided by Helix Sementes®, Brazil. Hybrids were obtained from crosses between 128 inbred lines and were evaluated for grain yield (GY) and plant height (PH). Field trials were carried out using a randomized complete block design with two replicates each, allocated across five sites for GY and three for PH during the growing season of 2015.

Inbred lines were genotyped via the Affymetrix® Axiom® Maize Genotyping Array (Unterseer et al. 2014) with 660K SNP markers. Quality control for SNPs was made based on call rate (CR), in which all markers with any missing data point were excluded, and minor allele frequency (MAF) procedures, in which markers with a low level of polymorphism (MAF < 0.05) were removed. Hybrid genotypes were scored by an allelic combination of homozygous markers of parental lines. After quality control, 37,625 SNP were used to compare the imputation methods.

Check the [CIMMYT DATA VERSE Repository](https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10887) for the full genotypic and phenotypic data

--------------------------------------------------------------------------------------------------------------------------------
<div id="p2" />

# USP data set

906 maize single-crosses obtained from a full dial- lel, according to Griffing’s method 4, divided into two heterotic groups, flint and dent, with 34 and 15 lines, respec- tively. Moreover, each heterotic group has a representative line, frequently used as the tester in our breeding program.

The experimental scheme used to evaluate the hybrids was an augmented block design (unreplicated trial) consisted of small blocks, each with 16 unique hybrids and two checks. Trials were carried out in Anhembi (22°50′51′′S, 48°01′06′′W, 466 m) and Piracicaba, at São Paulo State, Brazil (22°42′23′′S, 47°38′14′′W, 535 m), during the second growing season of 2016 and 2017, cultivated between January to June. In both sites and years, the hybrids were evaluated under two nitrogen (N) levels, low (LN) with 30 kg N ha−1, and normal (NN) with 100 kg N ha−1.

The genotyping of the 49 tropical inbred lines was per- formed by Affymetrix® platform, containing about 614,000 SNPs (Unterseer et al. 2014). Then, markers with low call rate (< 95%), minor allele frequency (MAF < 0.05) and heterozygous loci on at least one individual were removed. The missing markers were imputed using the [snpReady](https://github.com/italo-granato/snpReady) R package. Finally, the resulting 146,365 SNPs high-quality polymorphic SNPs were used to build the artificial hybrids genomic matrix, deduced by combining the genotypes from its two parents.

Check the [Mendeley Repository](https://data.mendeley.com/datasets/tpcw383fkm/3) for the full genotypic and phenotypic data


--------------------------------------------------------------------------------------------------------------------------------
<div id="p3" />

# Environmental Typing

The environmental typing (envirotyping) pipeline were conducted using the functions of the [EnvRtype](https://github.com/allogamous/EnvRtype) R package. This package has functions for supporting the collection of environmental data **get_weather()**, processing environmental data **processWTH()** and build of the **W** matrix of envirotype covariables **W.matrix()**. Finally, this package helps the construction of the genomic x envirotyping kinships using **get_kernels()**  functions, which easily can run into a Bayesian Genomic-enabled Prediction framework implemented in [BGGE](https://github.com/italo-granato/BGGE) package. Bellow we present a brief example of the use of both packages to run a genomic prediction considering reaction norms.

 ----------------------------------------------------------------------------------------------------------------
<div id="p9" />

# Authorship and Acknowledgments

For any questions, suggestions or corrections, please contact the author of this tutorial at **Germano Costa neto** germano.cneto@usp.br

We thank prof. Jaime Cuevas for his availability for Deep Kernel codes. For more details on this subject, contact him or look for more information in Cuevas et al. (2019).

We thank CIMMYT and the University of São Paulo (USP) for their support in the development of the project.

The theories exposed on this page are presented in the work of Costa Neto et al (not published yet) and were developed under the authorship of Germano Costa Neto, Jose Crossa and Roberto Fritsche Neto.
    


