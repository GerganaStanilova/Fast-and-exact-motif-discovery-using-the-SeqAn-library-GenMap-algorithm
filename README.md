<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->


<!-- PROJECT LOGO -->
<br />
<div align="center">
  
<h3 align="center">Fast and exact motif discovery <br />using the SeqAn library GenMap algorithm</h3>
  <p align="center">
    Finding transcription factor binding sites in protomer regions of gene families.
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
<br />
Motif discovery methods are of great significance because they can point to
meaningful biological conclusions. One application of motif discovery is identi-
fying DNA regulatory elements. In particular, transcription factor binding sites
since they are the target of drug therapy given their crucial role in modulating
transcription mechanisms.

The planted (l,d)-motif search is commonly used as a method for benchmark-
ing the performance of motif finding algorithms. It hunts for a l-character long
pattern with d mutations in a set of sequences with a fixed length.

This project compares the performance of two algorithms that are given different (l,d)-problems.

One of them uses the <a href="https://github.com/cpockrandt/genmap" target="_blank" rel="noopener noreferrer">GenMap algorithm</a>, from the <a href="https://docs.seqan.de/seqan/3-master-dev/index.html" target="_blank" rel="noopener noreferrer">SeqAn Library</a> and the other one is the implementation in C++ of the projection-algorithm from the paper 
 <a href="https://pubmed.ncbi.nlm.nih.gov/12015879/" target="_blank" rel="noopener noreferrer">"Finding motifs using random projections"</a> written by J. Buhler and M. Tompa.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



### Built With

<img height="50" src="https://docs.seqan.de/seqan/3-master-dev/seqan_logo.svg" /> &nbsp;&nbsp;&nbsp;&nbsp;<img height="50" src="https://upload.wikimedia.org/wikipedia/commons/1/18/ISO_C%2B%2B_Logo.svg" /> &nbsp;&nbsp;&nbsp;<img height="50" src="https://www.vectorlogo.zone/logos/cmake/cmake-ar21.svg" />


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

Install SeqAn3 following the instructions on the SeqAn <a href="https://docs.seqan.de/seqan/3-master-dev/usergroup0.html" target="_blank" rel="noopener noreferrer">Setup page</a>.

Install GenMap following the instructions on the GenMap <a href="https://github.com/cpockrandt/genmap#installation" target="_blank" rel="noopener noreferrer">Github Page</a>.



### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm.git
   ```
2. Set the CPP Compiler
	* for the GenMap Motiffinder
   ```sh
   /path/to/genmap_motiffinder/build cmake -DCMAKE_BUILD_TYPE=Release ../src
   ```	
	* for the Projection Motiffinder
   ```sh
   /path/to/projection_motiffinder/build cmake -DCMAKE_BUILD_TYPE=Release ../src
   ```
   * for creating synthethic datasets
   ```sh
   /path/to/create_datasets/build cmake -DCMAKE_BUILD_TYPE=Release ../src
   ```
   * for implanting a motif in random sequences
   ```sh
   /path/to/implant_motif/build cmake -DCMAKE_BUILD_TYPE=Release ../src
   ```
   
4. Compile the file
   * for the GenMap Motiffinder
   ```sh
   /path/to/genmap_motiffinder/build make
   ```	
	* for the Projection Motiffinder
   ```sh
   /path/to/projection_motiffinder/build make
   ```
   * for creating synthethic datasets
   ```sh
   /path/to/create_datasets/build make
   ```
	* for implanting a motif in random sequences
   ```sh
   /path/to/implant_motif/build make
   ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

1. The synthethic datasets in the "datasets" folder for each motiffinder can be used to run the program


* To run the program with synthetic data, execute:

   ```sh
   /path/to/build ./motiffinder ../path/to/the/synthetic/fasta/files  ../path/to/the/parameters/csv/file  ../path/to/the/correct/results/csv/files  numberofdatasets
   ```
   
* For example to run the program with 10 datasets in which there are hidden motifs with length 10 containing mutations with up to 2 mismatches, execute:

    ```sh
   /path/to/build ./motiffinder ../datasets/10-2/syn_synthetic_data_10_2.fasta ../datasets/10-2/parameters_10_2.csv ../datasets/10-2/syn_planted_motif_10_2.csv 10
   ```

	This screenshot shows the last part of the output after executing the last command:
	
	<img height="450" src="https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/blob/master/images/genmap_example.png" />


2. Alternatively, new synthetic datasets can be created by executing:


	```sh
   /path/to/create_datasets/build ./createdatasets t n l d
   ```
	where t is the number of the sequences, n is the length of one sequence, l is the length of the motif, and d is the maximum number of mutations the motif can have.
	
	For example, to create a dataset with 20 sequences of length 600 with a 10-character long motif with up to 2 mismatches, execute:

	```sh
   /path/to/create_datasets/build ./createdatasets 20 600 10 2
   ```
   You will get a fasta-file with the sequences and a csv-file with the correct plantet motif and its positions in the sequences located in the "build"-folder:
   
   <img height="50" src="https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/blob/master/images/syn_planted_motif_example.png" />


3. To run the program with biological data, execute:

   ```sh
   /path/to/build ./motiffinder ../path/to/the/synthetic/fasta/file  ../path/to/the/parameters/csv/file 
   ```
   
* In the folder "datasets/gcn4_biological" there is a file (GCN4_promoter_regions.fasta) that contains real experimentally verified promoter regions.

	For example running the program with the GCN4 promoter regions would look like this:

	```sh
   /path/to/build ./motiffinder ../datasets/gcn4_biological/GCN4_promoter_regions.fasta ../datasets/gcn4_biological/parameters_10_2.csv
   ```
   Since the program cannot account for the nucleotide background distribution, the results are suboptimal.
   
* To generate new random background sequences and implant the motif of GCN4, execute:

	```sh
   /path/to/implant_motif/build ./implantmotif t n 
   ```
   
   where t is the number of the sequences and n is the length of one sequence.
   
   
   
   
_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>





<!-- CONTACT -->
## Contact

Gergana Stanilova - stanilova.gergana@gmail.com

Project Link: [https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm](https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [GitGub README Template](https://github.com/othneildrew/Best-README-Template)
* []()
* []()

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-url]: https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm.svg?style=for-the-badge
[forks-url]: https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/network/members
[stars-shield]: https://img.shields.io/github/stars/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm.svg?style=for-the-badge
[stars-url]: https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/stargazers
[issues-shield]: https://img.shields.io/github/issues/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm.svg?style=for-the-badge
[issues-url]: https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/issues
[license-shield]: https://img.shields.io/github/license/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm.svg?style=for-the-badge
[license-url]: https://github.com/GerganaStanilova/Fast-and-exact-motif-discovery-using-the-SeqAn-library-GenMap-algorithm/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
[seqan-logo]: https://docs.seqan.de/seqan/3-master-dev/seqan_logo.svg
[seqan-url]: https://docs.seqan.de/seqan/3-master-dev/index.html
[cpp-logo]: https://upload.wikimedia.org/wikipedia/commons/1/18/ISO_C%2B%2B_Logo.svg
