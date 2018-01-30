FROM grandforest-web-common
MAINTAINER Simon J. Larsen <simonhffh@gmail.com>

COPY . /srv/shiny-server

CMD ["/usr/bin/shiny-server.sh"]
