## Contents of scripts/

Most of the files in the scripts/ folder are used for building OpenZGY
inside a docker container for various architectures. You don't need
these and you don't need docker if you just want to build a single
version of the code that will run on your current architecture.

That being said, the docker files might help setting up the
prerequisites correctly. And using docker avoids polluting your
current environment with a lot of extra packages. The downside is that
if you are not already familiar with docker then this approach might
involve more effort. And the scripts here need to do more and are thus
somewhat more complex than building by hand.
