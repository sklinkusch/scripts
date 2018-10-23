# scripts

This is a small description of programs within this folder

## fahrinfo

programs to read the data from fahrinfo.vbb.de, manipulate and display them in a perl-Tk window
* **fahrinfo-dual.pl**: read arrivals and departure of a station/bus stop
* **fahrinfo-duo-sort.pl**: read departures or arrivals of two stations/bus stops
* **fahrinfo-duo-sort-filter.pl**: read departures or arrivals of two stations/bus stops and filter them
* **fahrinfo-ffilter-ddual.pl**: read departures or arrivals of a station/bus stop and filter them according to two different schemes
* **fahrinfo-ffilter-dual.pl**: read departures and arrivals of a station/bus stop and filter them according to two different schemes
* **fahrinfo-filter-dual.pl**: read departures and arrivals of a station/bus stop and filter them
* **fahrinfo-filter.pl**: read departures or arrivals of a station/bus stop and filter them
* **fahrinfo-fstrict-dual.pl**: read departures and arrivals of a station/bus stop (w/o neighboring stops)
* **fahrinfo-fstrict-duo.pl**: read departures or arrivals of two different stations/bus stops (w/o neighboring stops)
* **fahrinfo-fstrict-filter.pl**: read departures or arrivals of a station/bus stop (w/o neighboring stops) and filter them
* **fahrinfo-fstrict.pl**: read departures or arrivals of a station/bus stop (w/o neighboring stops)
* **fahrinfo.pl**: read departures or arrivals of a station/bus stop
* **fahrinfo-strict.pl**: read departures or arrivals of a station/bus stop (with information on neighboring stops/track)

further programs
* **fahrinfo-dlinks.sh**: print timetable of fahrinfo.vbb.de into the command line (usually in less)
* **fahrinfo-elinks.sh**: show timetable of fahrinfo.vbb.de in elinks
* **fahrinfo-elinksd.sh**: print timetable of fahrinfo.vbb.de into the command line (not in less)
* **fahrinfo-w3m.sh**: show timetable of fahrinfo.vbb.de in w3m

data file
* **Fahrinfo.param**: including the update interval in seconds

routines file
* **Fahrinfo.pm**: routines for the fahrinfo.pl program suite
