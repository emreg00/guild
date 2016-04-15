latex2html -split 0 -show_section_numbers -local_icons -no_navigation -no_footnode -auto_prefix -no_subdir -info 0 guild.tex
vi -e guild.html < changes.vim

#scp guild.html ben-yehuda:/usr/local/apache2/htdocs/web/data/guild.software.html
#scp guild.pdf ben-yehuda:/usr/local/apache2/htdocs/data/guild/guild_handbook.pdf

