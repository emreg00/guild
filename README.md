# GUILD - Genes Underlying Inheritance Linked Disorders 

GUILD (Genes Underlying Inheritance Linked Disorders) is a framework built for the prioritization of disease candidate genes using a priori gene-disease associations and protein interactions. GUILD consists of implementations of 8 algorithms: NetScore, NetZcore, NetShort, NetCombo, fFlow, NetRank, NetWalk and NetProp. NetScore, NetZcore, NetShort, fFlow and NetRank are implemented in C++ while NetWalk and NetProp are implemented in R and NetCombo is a Python script combining the results of NetScore, NetZcore and NetShort. In this manual, we describe how to use these programs included in GUILD framework. 

In brief, GUILD is a software suite providing command line interface for the following network-based disease-gene prioritization algorithms:

* NetScore
* NetZcore
* NetShort
* NetCombo
* Functional Flow
* NetRank (PageRank with priors)
* Random walk with restart
* Network Propagation 

Documentation available at [doc](doc) and on the [web](http://sbi.imim.es/GUILD.php).

Contents:
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [Input files](#input-files)
* [Citation](#citation)

===

## Requirements
* GCC (GNU project C/C++ compiler) (version 4.3 or higher)
* make (GNU project utility to maintain groups of programs)
* BOOST lib (tested with boost_1_41_0 and boost_1_62_0)
* R (version 2.12.1 or higher) (only required for running NetWalk and NetProp)

Unix-like operating systems typically ship with these programs. If not these programs are freely available online. Note that, Windows users can have the fundamental environment for the installation (GCC and make) through MinGW (http:// http://www.mingw.org ) or Cygwin ( http://www.cygwin.com ). R is a free software environment for statistical computing and graphics and available at http://www.r-project.com /. 

For Python scripts used to create input files see guild_utilities in [toolbox](http://github.com/emreg00/toolbox) package.

## Installation

Clone this repository or download and unpack the source package [guild.tar.gz](http://sbi.imim.es/web/GUILD.php) e.g. as follows

Using github
```
\$> git clone https://github.com/emreg00/guild.git guild
```
If you use the source code from github, you have to download the [Boost C++ lib](http://www.boost.org/users/download/) and
extract it in "guild" folder where the repository is cloned.

Using SBI mirror 
```
\$> wget http://sbi.imim.es/data/guild/guild.tar.gz
\$> tar xvzf guild.tar.gz
```

Then, go to the cloned/extracted directory

```
\$> cd guild
```

Next, go to the src folder and issue make command as below. Beware that make command in MinGW can have a different name (e.g. mingw32-make.exe). 
**MacOS users:** use `make -f Makefile.mac` instead of `make` command below (available in the github version of the code).

```
\$> cd src 
\$> make
```

An executable named guild should be created under the "guild" folder. Try running it as follows

```
\$> cd ..
\$> ./guild
```

If you get the following output when you run it, the installation is successfully completed.

```
./guild [ Copyleft (GPLv3) - 2011 - Emre Guney (Universitat Pompeu Fabra) ]

 Arguments: 
     -s <prioritization_method>{NetScore:s|NetZcore:z|NetShort:d|fFlow:f|NetRank:r}
     -n <node_file>
     -e <edge_file>
     -o <output_file>
     -i <number_of_iterations>
     -r <number_of_repetitions>
     -t <seed_score_threshold>
     -x <number_of_sampled_graphs>
     -d <sampled_graph_prefix>
     -h
```

Otherwise make sure that you have recent versions of GCC and make installed, check the steps above and retry compiling. 

## Usage

### General overview

For algorithms implemented in C++, a typical GUILD call consist of several mandatory arguments (such as name of the input/output files and type of the prioritization method) followed by prioritization method specific arguments. Mandatory arguments common to all prioritization methods are explained below, method specific arguments are described in the later sections for each method separately. Possible arguments for a GUILD executable call is as follows:

```
\$> ./guild -s <prioritization_method> -n <node_file> -e <edge_file> -o <output_file> 
	    -i <number_of_iterations> -r <number_of_repetitions> -t <seed_score_threshold> 
	    -x <number_of_sampled_graphs> -d <sampled_graph_prefix>
```

where;

prioritization_method:
    The type of the prioritization algorithm, available values are

```
        s: NetScore
        z: NetZcore
        d: NetShort
        f: fFlow
        r: NetRank
```

For algorithms implemented in R (NetWalk and NetProp), a typical GUILD call would look like:

```
\$> R --slave --args <node_file> <edge_file> <output_file> <use_propagation> < random_walk.r
```

where all arguments are as explained before except "use_propagation", which -if provided- converts NetWalk algorithm to NetProp. 


### Method description and example runs 

* NetScore (suggested n_repetition: 3, n_iteration: 2)

NetScore adopts a message-passing scheme such that each node sends the information associated with it as a message to all of its neighbours and the neighbours convey these messages to their neighbours. NetScore takes into consideration alternative shortest paths within the distance of at most number_of_iterations links at each so called repetition-cycle. At the end of the repetition cycle, the node scores are updated according to messages recieved so far and the message passing is restarted. 

```
\$> ./guild -s s -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt -r 3 -i 2
```

* NetZcore (suggested n_iteration: 5, n_random_graph*: 100)

NetZcore assigns a normalized score using the distribution of the scores of neighbouring nodes. The normalization uses a random model of networks and it is calculated with the Z-score formulae: z=(x-m)/s, where m is the average of scores of neighbouring nodes with similar distribution in the random network and s is the standard deviation. The distribution is obtained with hundred network-replicates obtained by randomly shuffling the scores among nodes with similar degree (i.e. 100 random networks preserving the original topology). 

```
\$> ./guild -s z -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt -i 5 
	    -d data/test_interactions.txt. -x 100
```

An integer from 1 to n_random_graph will be appended at the end of prefix_of_random_graph (given by -d option) while reading random graphs (thus prefix_of_random_graphs should also include the directory under which graph files reside). A python script named "create_random_networks_for_netzcore.py" is provided for creating random networks that are going to be used by NetZcore. It requires Python (version 2.5.2 or higher) and Python NetworkX (version 1.1 or higher) package to be installed in your system. The following command would create 100 random networks with the same topology of given input network "data/test_interactions.txt" with the prefix of "data/test_interactions.txt." (appends a dot at the end of the provided egde scores file name).

```
\$> python src/create_random_networks_for_netzcore.py data/test_interactions.txt 100
```

* NetShort

NetShort accumulates the weighted shortest path lengths between a node and the rest of nodes in the network, where each edge-weight is inversely proportional to the average of the scores of the two nodes connected by the edge (i.e.edges connecting high scoring nodes are shorter). 
There is no method specific parameter for NetShort, however note that algorithm uses the phenotypic association scores in the edge scores file (edge_file) rather than the node scores file (e.g. the average of the scores of the nodes the edge in concern connects). A python script to create netshort specific edge scores file is provided for convenience (see below). 

```
\$> ./guild -s d -n data/test_proteins.txt -e data/test_interactions_for_netshort.txt -o output.txt
```

Edge scores are average relevance scores of the nodes they connect. A python script named "convert_network_for_netshort.py" is provided for creating edge scores file that is going to be used by NetShort (where original edge scores are multiplied by average of the scores of the nodes the edges belong to). It requires Python (version 2.5.2 or higher). The following command would convert the original edge scores file "data/test_interactions.txt" to a NetShort specific "data/test_interactions_for_netshort.txt" egde scores file using node scores information in "data/test_proteins.txt".

```
\$> python src/convert_network_for_netshort.py data/test_proteins.txt data/test_interactions.txt 
	data/test_interactions_for_netshort.txt
```

* Function Flow (suggested n_iteration: 5), reimplementation of [Nabieva et al.](http://www.ncbi.nlm.nih.gov/pubmed/15961472)

In fFlow (based on the algorithm of Functional Flow in Nabieva et al. 2005), at each iteration, annotation scores flow from nodes with higher score towards nodes with lower scores at the amount of the capacity of the edge through which the nodes are connected. 

```
\$> ./guild -s f -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt -i 5 -t 1.0
```

In this example all nodes that have a score equal or greater than 1.0 are considered seeds.

* NetRank (Page rank with priors), reimplemantation of [Chen et al.](http://www.ncbi.nlm.nih.gov/pubmed/19465376)

NetRank (based on the ToppGene algorithm proposed by Chen et al. 2009) uses Page Rank with priors algorithm (where a random surfer is more likely to end up in initially relevant nodes) to score a node in terms of phenotypic relevance. The damping factor is 0.15 and the number iterations for convergence is defined by the number_of_iterations parameter (see below). 
```

\$> ./guild -s r -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt
```

* NetWalk (Random walk with restart), reimplemantation of [Kohler et al.](http://www.ncbi.nlm.nih.gov/pubmed/18371930)

NetWalk (based on the Random walk with restarts algorithm proposed by Kohler, et al., 2008) iteratively simulates random transitions of a walker from a node to a randomly selected neighbour node and where at any time step the walk can be restarted depending on a predefined probability. Random walk with restarts is slightly different than PageRank with priors in the way that it normalizes the link weights. The convergence is decided by either having a probability difference less than 10e-6 between two consecutive time steps or achieving the limit of the number of iterations, set to 50 (though in practice less than 20 iterations are typically sufficient to satisfy the first criterion). There is no method specific parameter for NetWalk. 
The [random_walk.r](src/random_walk.r) is used to run NetWalk as follows:

```
\$> R --slave --args data/test_proteins.txt data/test_interactions.txt output.txt < random_walk.r
```

* NetProp (Network propagation), reimplementation of [Vanunu et al.](http://www.ncbi.nlm.nih.gov/pubmed/20090828)

NetProp (based on the Network propagation algorithm proposed by Vanunu, et al., 2010) modifies random walk with restarts such that the link weight is normalized not only by number of outgoing edges but also by number of incoming edges. The convergence is decided as it is done for NetWalk. There is no method specific parameter for NetProp. 
Similar to NetWalk, [random_walk.r](src/random_walk.r) is used with propagation=True argument, given as an integer in the command line.

```
\$> R --slave --args data/test_proteins.txt data/test_interactions.txt output.txt 1 < random_walk.r
```

The "1" turns "propagation" mode on.

* NetCombo

NetCombo combines the output scores from NetScore, NetZcore and NetShort in a consensus scheme by averaging normalized scores (z-scores) of a node in all of these methods. It requires the output files of NetScore, NetZcore and NetShort.
Therefore an example run for NetCombo is as follows: 

```
\$> python src/combine_scores.py output_netscore.txt output_netzcore.txt output_netshort.txt output.txt
```

## Input files

In a nutshell, the input files contain:
* node_scores: node id whitespace(s) node score
* edge_scores: node id whitespace(s) interaction score whitespace(s) node id
* edge_scores_as_node_scores: similar to edge_scores file but edge score is assigned averaging scores of the nodes of the edge in node_scores file

node_file:
    Input node scores file containing node (e.g. protein or gene) identifier followed by its phenotypic relevance score (e.g. association with the disease phenotype for that protein/gene) on each line. The values need to be separated by whitespace(s). That is;

```
    <node_id> <node_score>
```

edge_file:
    Input edge scores file containing node (e.g. protein or gene) identifier followed by score of the edge its phenotypic relevance score (e.g. association with the disease phenotype for the proteins/genes it is connecting) and node identifier (the interaction partner) on each line. The values are separated by whitespace(s). Thus, a line in this file looks like;

```
    <node_id> <edge_score> <node_id>
```

output_file:
    Output node scores file containing node (e.g. protein or gene) identifier followed by its calculated phenotypic relevance score (e.g. association with the disease phenotype for that protein/gene) on each line. The values are separated by whitespace(s). The format of a line would be;

```
    <node_id> <node_score>
```

## Citation

Guney E, Oliva B. Exploiting Protein-Protein Interaction Networks for Genome-Wide Disease-Gene Prioritization. PLoS ONE 7(9): e43557 (2012). [link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043557)

