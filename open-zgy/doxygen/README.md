This folder contains a snapshot of the OpenZGY documentation, but only as pdf files and only for the public api. It might not be completely up to date.

The latest version of the OpenZGY documentation is built by Doxygen and is available after the project has been built. If you are building the software yourself then you will find the documentation in build/deploy/native.

If you don't want to build yourself then you can fetch the documentation from the CI/CD pipeline in gitlab as follows.

* Open the [OpenZGY project page in OSDU](https://community.opengroup.org/osdu/platform/domain-data-mgmt-services/seismic/open-zgy)
* Go to CI/CD -> Pipelines
* Select the latest good build on master and click on the "passed" box.
* Click on the doc-compile box
* Click "Browse"
* Navigate to build/deploy/native or build/deploy/pure
* Download the documentation you want.

apidoc has the public api only. intdoc also includes implementation.

The pdf version is good for a quick look and you can read it directly on the web. The html version is much better organized and has good search capabilities but you need to download and unpack it before viewing it locally. The OSDU web server cannot serve these pages directly. It can only serve single page documents.

At some point the documentation built on gitlab might be deployed to a proper web server that allows the html version to be displayed online.
