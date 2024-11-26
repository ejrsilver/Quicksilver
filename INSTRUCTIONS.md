# Instructions for Configuring and Running Quicksilver

Set the following environment variable with the path of your PowerAPI install.
This is what it looks like on my account on the CAC machine:

```bash
POWER_LOC=/global/home/hpc5574/pwrapi-ref/build/install
```

Also edit the [xml_profile](./xml_profile) with the path to the xml file you want to use.

Then, run `source xml_profile`. This will configure the environment variables.

To build Quicksilver, `cd` into [src](./src/) and build with `make`.

To run the executable, `cd` back to the base directory and run `./src/qs`.
