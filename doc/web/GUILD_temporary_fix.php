<?php
        $root_path = dirname(__FILE__);
	
	function curPageURL() {
		$request_uri = $_SERVER["REQUEST_URI"];
		$temp = explode("?",$request_uri);
		$request_uri = $temp[0];
	
 		$pageURL = 'http';
		 if ($_SERVER["HTTPS"] == "on") {$pageURL .= "s";}
			 $pageURL .= "://";
		 if ($_SERVER["SERVER_PORT"] != "80") {
			 $pageURL .= $_SERVER["SERVER_NAME"].":".$_SERVER["SERVER_PORT"].$request_uri;
		 } else {
			$pageURL .= $_SERVER["SERVER_NAME"].$request_uri;
		 }
		 return $pageURL;
	}

	$current_page_url = curPageURL();
	if( isset($_GET["page"]) ){
	        $pagenum=$_GET["page"];
	}
	else{
		$pagenum="guild";
	}
?>
<html>
<head>
<title> GUILD (Genes Underlying Inheritance Linked Disorders) </title>
</head>
<body>

<h1>GUILD (Genes Underlying Inheritance Linked Disorders)</h1>

<p><span>The past decade has witnessed dramatic advances in genome sequencing and a substantial shift in the number of genome wide association studies (GWAS). These efforts have expanded considerably our knowledge on the sequential variations in Human DNA and their consequences on the human biology. Complex genetic disorders often involve products of multiple genes acting cooperatively. Nevertheless, pinpointing the decisive elements of such disease pathways is still a challenge. Recently, network biology has proven its use in identifying candidate genes associated with a disease based on the simple observation that proteins translated by phenotypically related genes tend to interact, so called guilt-by-association principle. Here, we present GUILD (Genes Underlying Inheritance Linked Disorders), a network-based prioritization framework to unveil genes associated with a disease phenotype (disease-genes). In GUILD, we exploit several communication mechanisms between disease-genes emerging from the topology of the interaction network. We used three sources of gene-phenotypic association to specify nodes involved in a disorder (including Online Mendelian Inheritance in Man database)previously published data sets. Analyses on multiple human disease phenotypes demonstrated that the methods proposed in GUILD effectively prioritize genes even when the linkage information is not known <em> a priori</em>. We compared the algorithms we have developed with state-of-the-art prioritization methods such as PageRank with priors, Functional-Flow, Random walk with restart and Network propagation. We also tested the robustness of the approaches proving the effect of the network properties and the independence with the number of original genes/proteins associated with the function or phenotype. Finally, we applied GUILD to prioritize genes in the case of Alzheimer Disease.   <!--
<div align="center"><img src="galeries/guild/overview.png" _mce_src="galeries/guild/overview.png" alt="GUILD overview" /></div>
--> </span></p>
<p>GUILD (Genes Underlying Inheritance Linked Disorders) is a framework built for the prioritization of disease candidate genes using a priori gene-disease associations and protein interactions. GUILD consists of implementations of 7 algorithms: NetScore, NetZcore, NetShort, fFlow, NetRank, NetWalk and NetProp. NetScore, NetZcore, NetShort, fFlow and NetRank are implemented in <em>C++</em> while NetWalk and NetProp are implemented in <em>R</em>. In this manual, we describe how to use these programs included in GUILD framework.</p>
<p> </p>
<h3>Download</h3>
<p>The software used in GUILD framework is distributed freely under GNU Public License. Download <a href="data/guild/guild.tar.gz"> the source code, example input files and the software manual here</a>. See handbook either in the downloaded achieve (under "doc" directory) or below for the instructions on installation and usage. A <a href="data/guild/guild_handbook.pdf"> PDF version of the handbook</a> is also provided for convenience.</p>
<h3>Handbook</h3>
<h4><a name="SECTION00020000000000000000"> 1 Requirements</a></h4>
<p> </p>
<ul>
<li>GCC (GNU project C/C++ compiler) (version 4.3 or higher) </li>
<li>make (GNU project utility to maintain groups of programs) </li>
<li>R (version 2.12.1 or higher) <em>(only required for running NetWalk and NetProp)</em> </li>
</ul>
<p>Unix-like operating systems typically ship with these programs. If not these programs are freely available online. Note that, Windows users can have the fundamental environment for the installation (<em>GCC</em> and <em>make</em>) through <em>MinGW</em> (http://<a href="http://www.mingw.org"> http://www.mingw.org </a>) or <em>Cygwin</em> (<a href="http://www.cygwin.com"> http://www.cygwin.com </a>). <em>R</em> is a free software environment for statistical computing and graphics and available at <a href="http://www.r-project.org"> http://www.r-project.com </a>/.</p>
<p> </p>
<h4><a name="SECTION00030000000000000000"> 2 Installation</a></h4>
<p>Download and unpack the source package <em>guild.tar.gz</em> located at <a href="GUILD.php"> http://sbi.imim.es/GUILD.php </a> e.g. as follows</p>
<p> </p>
<pre>{bash}
\$&gt; tar xvzf guild.tar.gz
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> Then, go to the extracted directory</p>
<pre>{bash}
\$&gt; cd guild
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> Next, go to the src folder and issue <em>make</em> command as below. Beware that <em>make</em> command in MinGW can have a different name (e.g. <em>mingw32-make.exe</em>)</p>
<pre>{bash}
\$&gt; cd src
\$&gt; make
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> An executable named <em>guild</em> should be created under the ``guild'' folder. Try running it as follows</p>
<p> </p>
<pre>{bash}
\$&gt; cd ..
\$&gt; ./guild
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> If you get the following output when you run it, the installation is successfully completed.</p>
<p> </p>
<pre>{bash}
./guild [ Copyleft (GPLv3) - 2011 - Emre Guney (Universitat Pompeu Fabra) ]

 Arguments:
     -s &lt;prioritization_method&gt;{NetScore:s|NetZcore:z|NetShort:d|fFlow:f|NetRank:r}
     -n &lt;node_file&gt;
     -e &lt;edge_file&gt;
     -o &lt;output_file&gt;
     -i &lt;number_of_iterations&gt;
     -r &lt;number_of_repetitions&gt;
     -t &lt;seed_score_threshold&gt;
     -x &lt;number_of_sampled_graphs&gt;
     -d &lt;sampled_graph_prefix&gt;
     -h
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> Otherwise make sure that you have recent versions of <em>GCC</em> and <em>make</em> installed, check the steps above and retry compiling.</p>
<p> </p>
<h4><a name="SECTION00040000000000000000"> 3 Usage</a></h4>
<p>For algorithms implemented in C++, a typical GUILD call consist of several mandatory arguments (such as name of the input/output files and type of the prioritization method) followed by prioritization method specific arguments. Mandatory arguments common to all prioritization methods are explained below, method specific arguments are described in the later sections for each method separately. Possible arguments for a GUILD executable call is as follows:</p>
<p> </p>
<pre>{bash}
\$&gt; ./guild -s &lt;prioritization_method&gt; -n &lt;node_file&gt; -e &lt;edge_file&gt; -o &lt;output_file&gt;
            -i &lt;number_of_iterations&gt; -r &lt;number_of_repetitions&gt; -t &lt;seed_score_threshold&gt;
            -x &lt;number_of_sampled_graphs&gt; -d &lt;sampled_graph_prefix&gt;
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> where;</p>
<p> </p>
<dl> <dt><strong>prioritization_method:</strong></dt> <dd>The type of the prioritization algorithm,         available values are     
<ul>
<li>s: NetScore </li>
<li>z: NetZcore </li>
<li>d: NetShort </li>
<li>f: fFlow </li>
<li>r: NetRank </li>
</ul>
<p> </p>
</dd> <dt><strong>node_file:</strong></dt> <dd>Input node scores file containing node (e.g. protein or         gene) identifier followed by its phenotypic relevance score (e.g.         association with the disease phenotype for that protein/gene) on each         line. The values need to be separated by whitespace(s). That is;
<pre>{text}
&lt;node_id&gt; &lt;node_score&gt;
</pre>
<span style="color: #bfbfbf;"><span>text</span></span>
<p> </p>
</dd> <dt><strong>edge_file:</strong></dt> <dd>Input edge scores file containing node (e.g. protein or         gene) identifier followed by score of the edge its phenotypic         relevance score (e.g. association with the disease phenotype for the         proteins/genes it is connecting) and node identifier (the interaction         partner) on each line. The values are separated by whitespace(s). Thus,         a line in this file looks like;
<pre>{text}
&lt;node_id&gt; &lt;edge_score&gt; &lt;node_id&gt;
</pre>
<span style="color: #bfbfbf;"><span>text</span></span>
<p> </p>
</dd> <dt><strong>output_file:</strong></dt> <dd>Output node scores file containing node (e.g. protein         or gene) identifier followed by its ``calculated'' phenotypic         relevance score (e.g. association with the disease phenotype for that         protein/gene) on each line. The values are separated by whitespace(s).         The format of a line would be;
<pre>{text}
&lt;node_id&gt; &lt;node_score&gt;
</pre>
<span style="color: #bfbfbf;"><span>text</span></span>
<p> </p>
</dd> </dl>
<p>For algorithms implemented in R (NetWalk and NetProp), a typical GUILD call would look like:</p>
<p> </p>
<pre>{bash}
\$&gt; R --slave --args &lt;node_file&gt; &lt;edge_file&gt; &lt;output_file&gt; &lt;use_propagation&gt; &lt; random_walk.r
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> where all arguments are as explained before except ``use_propagation'', which -if provided- converts NetWalk algorithm to NetProp.</p>
<p> </p>
<h4><a name="SECTION00041000000000000000"> 3.1 NetScore</a></h4>
<p>NetScore adopts a message-passing scheme such that each node sends the information associated with it as a message to all of its neighbours and the neighbours convey these messages to their neighbours. NetScore takes into consideration alternative shortest paths within the distance of at most number_of_iterations links at each so called repetition-cycle. At the end of the repetition cycle, the node scores are updated according to messages recieved so far and the message passing is restarted.</p>
<p><img src="galeries/guild/guild-img1.png" border="0" alt="includegraphics[scale=0.5]{netscore.eps}" width="384" height="215" align="BOTTOM" /></p>
<p><br /> <br /> Method specific parameters for NetScore are;</p>
<dl> <dt><strong>number_of_repetitions:</strong></dt> <dd>The number of resets (updating the scores         of the nodes to calculated scores so far) in the algorithm. Defines the         reach of the method (number of links to look further while calculating         score) in accordance with the number_of_iteration parameter. </dd> <dt><strong>number_of_iterations:</strong></dt> <dd>Number of iterations to apply the scoring         step of the prioritization algorithm. That is the length of shortest         paths to consider from a node to other nodes (messages of the nodes         that are this many links further are considered). </dd> </dl>
<p>The following is an example call to run NetScore algorithm using node file <em>node_score.txt</em> and edge <em>file edge_score.txt</em> with number of iteration and repetition parameters of <em>2</em> and <em>3</em> respectively, writing the calculated scores to a file named output.txt.</p>
<p> </p>
<pre>{bash}
\$&gt; ./guild -s s -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt -r 3 -i 2
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p> </p>
<h4><a name="SECTION00042000000000000000"> 3.2 NetZcore</a></h4>
<p>NetZcore assigns a normalized score using the distribution of the scores of neighbouring nodes. The normalization uses a random model of networks and it is calculated with the Z-score formulae: z=(x-m)/s, where m is the average of scores of neighbouring nodes with similar distribution in the random network and s is the standard deviation. The distribution is obtained with hundred network-replicates obtained by randomly shuffling the scores among nodes with similar degree (i.e. 100 random networks preserving the original topology).</p>
<p><img src="galeries/guild/guild-img2.png" border="0" alt="includegraphics[scale=0.5]{netzcore.eps}" width="466" height="514" align="BOTTOM" /></p>
<p><br /> <br /> Method specific parameters for NetZcore are;</p>
<dl> <dt><strong>number_of_iterations:</strong></dt> <dd>Number of iterations to apply the scoring         step of the prioritization algorithm. </dd> <dt><strong>number_of_sampled_graphs:</strong></dt> <dd>Number of the sampled networks (random         networks with similar characteristics -e.g. topology or degree         distribution- of original network) used for calculating expected mean         and standard deviation of scores by random. </dd> <dt><strong>sampled_graph_prefix:</strong></dt> <dd>The full prefix of the sampled networks. An         integer from 1 to n_sampled_graphs will be appended at the end of         this text while reading sampled networks (thus it should include the         directory under which random network files reside). Note that the         program itself does not create random networks and looks for already         existing random networks residing under the path given by this parameter.         A python script is provided for creating such networks (see below). </dd> </dl>
<p>An example call to run NetZcore algorithm where <em>data</em> folder contains 100 randomly generated networks starting with the prefix <em>test_interactions.txt.</em> (e.g. <em>test_interactions.txt.1</em>, <em>test_interactions.txt.2</em>, ..., <em>test_interactions.txt.100</em> is as follows:</p>
<p> </p>
<pre>{bash}
\$&gt; ./guild -s z -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt -i 5
            -d data/test_interactions.txt. -x 100
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> A python script named ``create_random_networks_for_netzcore.py'' is provided for creating random networks that are going to be used by NetZcore. It requires <em>Python</em> (version 2.5.2 or higher) and Python <em>NetworkX</em> (version 1.1 or higher) package to be installed in your system. The following command would create 100 random networks with the same topology of given input network ``data/test_interactions.txt'' with the prefix of ``data/test_interactions.txt.'' (appends a dot at the end of the provided egde scores file name).</p>
<p> </p>
<pre>{bash}
\$&gt; python src/create_random_networks_for_netzcore.py data/test_interactions.txt 100
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p> </p>
<h4><a name="SECTION00043000000000000000"> 3.3 NetShort</a></h4>
<p>NetShort accumulates the weighted shortest path lengths between a node and the rest of nodes in the network, where each edge-weight is inversely proportional to the average of the scores of the two nodes connected by the edge (i.e.edges connecting high scoring nodes are shorter).</p>
<p><img src="galeries/guild/guild-img3.png" border="0" alt="includegraphics[scale=0.5]{netshort.eps}" width="384" height="203" align="BOTTOM" /></p>
<p><br /> <br /> There is no method specific parameter for NetShort, however note that algorithm uses the phenotypic association scores in the edge scores file (<em>edge_file</em>) rather than the node scores file (e.g. the average of the scores of the nodes the edge in concern connects). A python script to create netshort specific edge scores file is provided for convenience (see below).</p>
<p><br /> <br /> Thus an example NetShort run would be;</p>
<p> </p>
<pre>{bash}
\$&gt; ./guild -s d -n data/test_proteins.txt -e data/test_interactions_for_netshort.txt -o output.txt
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p><br /> <br /> A python script named ``convert_network_for_netshort.py'' is provided for creating edge scores file that is going to be used by NetShort (where original edge scores are multiplied by average of the scores of the nodes the edges belong to). It requires <em>Python</em> (version 2.5.2 or higher). The following command would convert the original edge scores file ``data/test_interactions.txt'' to a NetShort specific ``data/test_interactions_for_netshort.txt'' egde scores file using node scores information in ``data/test_proteins.txt''.</p>
<p> </p>
<pre>{bash}
\$&gt; python src/convert_network_for_netshort.py data/test_proteins.txt data/test_interactions.txt
        data/test_interactions_for_netshort.txt
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p> </p>
<h4><a name="SECTION00044000000000000000"> 3.4 fFlow</a></h4>
<p>In fFlow (based on the algorithm of Functional Flow in Nabieva et al. 2005), at each iteration, annotation scores flow from nodes with higher score towards nodes with lower scores at the amount of the capacity of the edge through which the nodes are connected.</p>
<p><br /> <br /> The method specific parameters for fFlow are;</p>
<dl> <dt><strong>number_of_iterations:</strong></dt> <dd>Number of iterations to apply the scoring         step of the prioritization algorithm. </dd> <dt><strong>seed_score_threshold:</strong></dt> <dd>All the nodes that have higher than this         given threshold score will be considered as seeds for the method and         assigned infinite scores during scoring. </dd> </dl>
<p>An example call of fFlow where all nodes that have a score higher than <img src="galeries/guild/guild-img4.png" border="0" alt="$1.0$" width="27" height="18" align="BOTTOM" /> are seeds is;</p>
<p> </p>
<pre>{bash}
\$&gt; ./guild -s f -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt -i 5 -t 1.0
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p> </p>
<h4><a name="SECTION00045000000000000000"> 3.5 NetRank</a></h4>
<p>NetRank (based on the ToppGene algorithm proposed by Chen et al. 2009) uses Page Rank with priors algorithm (where a random surfer is more likely to end up in initially relevant nodes) to score a node in terms of phenotypic relevance. The damping factor is 0.15 and the number iterations for convergence is defined by the number_of_iterations parameter (see below).</p>
<p><br /> <br /> The method specific parameter for NetRank is;</p>
<dl> <dt><strong>number_of_iterations:</strong></dt> <dd>Number of iterations to apply the scoring step of the prioritization algorithm. In the case of Page Rank algorithm it defines the number of the times for calculating page ranks for the nodes (convergence criterion for the algorithm). This number is set to 20 by default. </dd> </dl>
<p>Therefore an example run for NetRank is as follows:</p>
<p> </p>
<pre>{bash}
\$&gt; ./guild -s r -n data/test_proteins.txt -e data/test_interactions.txt -o output.txt
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span></p>
<p> </p>
<h4><a name="SECTION00046000000000000000"> 3.6 NetWalk</a></h4>
<p>NetWalk (based on the Random walk with restarts algorithm proposed by Kohler, et al., 2008) iteratively simulates random transitions of a walker from a node to a randomly selected neighbour node and where at any time step the walk can be restarted depending on a predefined probability. Random walk with restarts is slightly different than PageRank with priors in the way that it normalizes the link weights. The convergence is decided by either having a probability difference less than 10e-6 between two consecutive time steps or achieving the limit of the number of iterations, set to 50 (though in practice less than 20 iterations are typically sufficient to satisfy the first criterion). There is no method specific parameter for NetWalk.</p>
<p><br /> <br /> Therefore an example run for NetWalk is as follows:</p>
<p> </p>
<h4><a name="SECTION00047000000000000000"> 3.7 NetProp</a></h4>
<p>NetProp (based on the Network propagation algorithm proposed by Vanunu, et al., 2010) modifies random walk with restarts such that the link weight is normalized not only by number of outgoing edges but also by number of incoming edges. The convergence is decided as it is done for NetWalk. There is no method specific parameter for NetProp.</p>
<p><br /> <br /> Therefore an example run for NetProp is as follows:</p>
<p> </p>
<pre>{bash}
\$&gt; R --slave --args data/test_proteins.txt data/test_interactions.txt output.txt 1 &lt; random_walk.r
</pre>
<p><span style="color: #bfbfbf;"><span>bash</span></span> </p>
<h2>Acknowledgements</h2>
<p>The authors are indebted to the members of the SBI laboratory for their involvement in the fruitful discussions on the project.</p>
<p>Emre Guney is grateful to the financial support from <strong>"Departament d'Educació i Universitats de la Generalitat de Catalunya i del Fons Social Europeu"</strong>.</p>
<p>This work was supported by grants from <strong>Spanish Ministry of Science and Innovation (MICINN) BIO2008-0205</strong>, and with <strong>FEDER support</strong>, <strong>PSE-0100000-2007</strong> and <strong>PSE-0100000-2009</strong>. </p>

</body>
</html>

