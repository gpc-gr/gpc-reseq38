#
#
#

NAME := registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-shapeit
TAG := bcftools_1.11-shapeit_4.4.2

.PHONY: build save clean

docker:
	docker build -t $(NAME):$(TAG) .

singularity:
	mkdir -p ./build/singularity
	docker run \
		-it \
		--rm \
		--privileged \
		--volume /var/run/docker.sock:/var/run/docker.sock \
		--volume $(shell pwd)/build/singularity:/output \
		quay.io/singularity/docker2singularity \
			$(NAME):$(TAG)

clean:
	rm -rf ./build || true
