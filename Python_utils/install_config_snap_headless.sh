wget http://step.esa.int/downloads/7.0/installers/esa-snap_sentinel_unix_7_0.sh
printf 'o\n2\n2\n/opt/snap\n2,3,4,5\ny\n/usr/local/bin\ny\n\ny' | bash ./esa-snap_sentinel_unix_7_0.sh
rm ./esa-snap_sentinel_unix_7_0.sh
snap -- jdkhome  /usr/lib/jvm/java-1.11.0-openjdk-amd64