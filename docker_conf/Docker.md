# ProPIP docker 


## Get ProPIP 

```
git clone https://github.com/acg-team/ProPIP.git

```

## Build `local/propip` docker image

```
cd ProPIP
docker build -t local/propip -f docker_conf/Dockerfile .
```

## Test build

```
docker run --rm local/propip ProPIP 
```

## Run test data using `local/propip` docker image.

- ```-v `pwd`/tests:/tests``` bind local tests folder to container's /test folder
- `-w /tests`set container's working dir to /test

```
docker run --rm -v `pwd`/tests:/tests -w /tests  local/propip ProPIP params=/tests/input/test_01/conf
docker run --rm -v `pwd`/tests:/tests -w /tests  local/propip ProPIP params=/tests/input/test_02/conf
docker run --rm -v `pwd`/tests:/tests -w /tests  local/propip ProPIP params=/tests/input/test_03/conf
docker run --rm -v `pwd`/tests:/tests -w /tests  local/propip ProPIP params=/tests/input/test_04/conf
docker run --rm -v `pwd`/tests:/tests -w /tests  local/propip ProPIP params=/tests/input/test_05/conf
docker run --rm -v `pwd`/tests:/tests -w /tests  local/propip ProPIP params=/tests/input/test_06/conf
```


