for i in *.php
do
	echo $i
	vi -e $i < modify_links.vim
	#mv $i ${i%.*}.php 
done
