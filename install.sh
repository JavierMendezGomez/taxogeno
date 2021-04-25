#!/bin/sh
sudo mkdir /opt/taxogeno/
sudo useradd taxogeno -U -d /opt/taxogeno
sudo chown -R taxogeno:taxogeno /opt/taxogeno
sudo cp -r ./* /opt/taxogeno

cd /opt/taxogeno
sudo -u postgres psql -c "CREATE USER taxogeno;" 
sudo -u postgres psql -c "CREATE DATABASE taxogeno WITH OWNER taxogeno;"

sudo -u postgres pg_restore -d taxogeno /opt/taxogeno/sql/taxogeno_full.sql
sudo -u postgres pg_restore -d taxogeno /opt/taxogeno/sql/biosql_pgsql.sql

sudo -u taxogeno Rscript /opt/taxogeno/taxogeno_init.R


