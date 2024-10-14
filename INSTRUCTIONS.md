# Instructions for Configuring and Running Quicksilver

First, change the path of your Power API install in the [Makefile](./src/Makefile) to your install directory.

Then, run `source xml_profile`. This will make all necessary environment variables available, and allow running Power API commands.

Using the Power API command `pwrgen`, a file `hwloc.xml` will be generated. This will reflect your system configuration, without many details.

Replace the empty `<Plugins/>` tag with the following:
```xml
    <Plugins>
        <plugin name="RAPL" lib="libpwr_rapldev"/>
    </Plugins>
```

Then `cd` into [src](./src/) and build with `make`.

Next, to run with RAPL, use the following command: `sudo --preserve-env ./src/qs`.
