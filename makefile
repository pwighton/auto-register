docker: docker/dockerfile.areg	
	docker build -f ./docker/dockerfile.areg -t areg ./docker

docker-bash: docker
	docker run -it --rm areg /bin/bash
