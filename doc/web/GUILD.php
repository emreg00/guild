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

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>GUILD - <?php print $pagenum; ?> - Structural Bioinformatics Laboratory</title>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />
<meta http-equiv='Content-Language' content='en_US' />
<meta http-equiv='Content-Script-Type' content='text/javascript' />
<meta http-equiv='Content-Style-Type' content='text/css' />
<meta name='keywords' content='Structural Bioinformatics GRIB IMIM UPF Biomedical Informatics' />
<meta name='description' content='Structural Bioinformatics Group in GRIB (Research Group on Biomedical Informatics' />
<meta name='author' content='Emre Guney' />
<meta name='Robots' content='index,follow' />
<meta http-equiv='imagetoolbar' content='no' /><!-- disable IE's image toolbar -->
<link rel='stylesheet' type='text/css' href='templates/modified_10news/style.css' />
<link rel='stylesheet' type='text/css' href='css/lightneasy.css' />
<link rel='stylesheet' type='text/css' href='data/biana.css' />
</head>

<body>

	<h1> <img WIDTH="20%" HEIGHT="20%" src="http://sbi.imim.es/guild/guild.png"> Genes Underlying Inheritance Linked Disorders </h1>

	     <ul>
		<li><a href="<?php echo $current_page_url;?>?page=guild.introduction">Introduction</a></li>
		<li><a href="<?php echo $current_page_url;?>?page=guild.software">Software</a></li>
		<li><a href="<?php echo $current_page_url;?>?page=guild.acknowledgements">Acknowledgements</a></li>
		<li><a href="<?php echo $current_page_url;?>?page=guild.citation">Citation</a></li>
	     </ul>
	<?php
		if( $pagenum=="guild" ){
			include("guild.introduction.php");
			print "<br/><hr />";
			include("guild.software.php");
			print "<br/><hr />";
			include("guild.acknowledgements.php");
			print "<br/><hr />";
			include("guild.citation.php");
		}
		else{
			include($pagenum.".php");	
		}

	?>
	<br />
</body>
</html>

