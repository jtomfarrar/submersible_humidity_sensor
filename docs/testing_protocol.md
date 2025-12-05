June 13, 2025

SPdeS comments Aug 18, 2025

# Submersible humidity sensor testing protocol

Major issues from previous tests:

(a) Probes not thermally equilibrated to water. They are not always
    submerged long enough to equilibrate.

(b) Not clear what true temperature and humidity actually is. No
    reference instrument was used to measure them.

(c) Probes did not reach same temperature during submerged time,
    suggesting (a) lack of equilibration, but we're not sure because (b)
    there was no reference measurement.

(d) Adjustment after exposure to air is confusing. This is exacerbated
    by (a); and (b) would provide clarity about the adjustment.

Link to [spreadsheet of previous
tests](https://docs.google.com/spreadsheets/d/1vdw-wr8KG5HfAMqWTgq3-DX9OUjbawgZ6YE9wdCr4aQ/edit?gid=1633130170#gid=1633130170)

**Dedicate a calibrated stand-alone ASIMET WXT to "live" with the SHS
probe and be used to provide references for SHS tests.**

## Revised protocol (06/13/2025):

1.  Use an ASIMET WXT as a reference standard for humidity and
    temperature (and, if outside, at least note wind, sunlight,
    weather).

2.  Sync WXT module clock and SHS clock to common standard

3.  Note time, turn on WXT

4.  Note time, turn on SHS logger, leave in air for \>5 minutes

5.  Note time, submerge for at least **2 minutes**

6.  Note time, remove from water for at least **30 minutes** (at least
    for next lab test)

    a.  SPdeS comment: 30 minutes should usually be long enough to
        equilibrate to wet bulb temperature, but is probably not long
        enough for the wet sensor to dry. We also need to characterize
        whether the wet-sock sensor dries off significantly, as this
        will affect interpretation of the two temperatures. We will be
        able to do this quicker and with more certainty with the
        reference measurements from the WXT.

7.  Repeat (5) and (6) two more times

## Data processing

-   Produce a netcdf file from the SHS with distinguishable variable
    names and unique name for the test (e.g. the date).

-   Include WXT data in the SHS netcdf file

    -   NetCDF files can't have 2 record dimensions, so 2 files might be
        better able to accommodate the two separate time dimensions for
        the SHS and the WXT.

## Preliminary tests with longer adjustment time periods

### Tom Farrar wrote:
It would be nice to exceed the "at least" times listed for the various transitions in this first run, especially the 2-min wet and 30-min dry times.  (This will let us understand what's really required in future runs.)  -Tom

### de Szoeke, Simon wrote:

It is a good idea to do some tests to determine the equilibration times, to inform the times in the protocol, which are actually guesses at periods that exceed equilibration times.

- We know the sensors are equilibrated when they change little, and fluctuate about a stationary mean, for a long time. We can gain confidence in this by looking at long time series.

- The equilibration time is one of the first things we need in order to interpret the tests. Once we know the time scale, we can refine the times in the protocol and do more standard (probably shorter) tests.

- The fluctuations after the temperature has equilibrated to the mean are useful, as they tell us something about the turbulence, independent of the sensor adjustment.

- It would also be useful to know how long the wet sensor stays wet, and what Twet does as the wick dries out. This could be done by running the last test overnight.

All these are benefits of doing some long tests. How long? We could try

5 min equilibration in water and 
2 hours in air

Ventilation will cut down on the time needed for sensor equilibration. I assume ambient ventilation "wind" will be measured by the WXT.
