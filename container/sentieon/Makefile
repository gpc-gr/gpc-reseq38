#
#
#

NAME := registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-sentieon
TAG := bcftools_1.11-sentieon_202010.02

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
