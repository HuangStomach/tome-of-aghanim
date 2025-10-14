docker run --detach --name mariadbA --env MARIADB_ROOT_PASSWORD=password --env MARIADB_DATABASE=testA --env MARIADB_USER=testA --env MARIADB_PASSWORD=balabala -p 3306:3306 mariadb:lts-noble
docker run --detach --name mariadbB --env MARIADB_ROOT_PASSWORD=password --env MARIADB_DATABASE=testB --env MARIADB_USER=testB --env MARIADB_PASSWORD=balabala -p 3307:3306 mariadb:lts-noble

