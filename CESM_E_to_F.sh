#####################################################################################
# A script for making boundary files for AMIP simulations with ACCESS ESM1.5
# Modification regarding using other SST as inputs are instructed
# Please make sure that the sea ice and SST have the same and standard resolution
# The reference is here:
#         http://www.cesm.ucar.edu/models/cesm1.2/cesm/doc/usersguide/x2306.html
#####################################################################################

#######################################################
# Setting Parameters
#######################################################

# SW: Gadi specific - load the required modules
module load intel-compiler/2021.8.0
module load netcdf
module load nco/5.0.5

# directory of this file and the bcgen folder
export dir_tool="/home/565/sw6175/prepping_forcings/CESM_BC"

# directory of SST and seaice data
export dir_data="/g/data/w40/sw6175/diddling/SST_patterns/Regridded/combined"

# Directory for diddling output
export dir_output="/g/data/w40/sw6175/diddling/SST_patterns/Diddled"

# Directory for temporary files
export dir_temp="/g/data/w40/sw6175/diddling/SST_patterns/Temp"

# Name of input file, asssumed to lie in dir_data
export input_file="regridded_input4mips_combined_floattime.nc"



# SPECIFY the first and last month and year for period in which observed 
# monthly mean data will be read. (The first month read must not preceed mon1, 
# iyr1, and the last month must not follow monn, iyrn).  

export mon1rd=1
export iyr1rd=1870

export monnrd=12
export iyrnrd=2016

# SPECIFY first and last month and year for entire period that will
# be treated (i.e., the period of interest plus buffers at ends).
# Note that the entire period treated should be an integral   
# number of years.  

# Try adding a buffer of one year on either side
export mon1=1
export iyr1=1869
export monn=12
export iyrn=2017

# SPECIFY the first and last month and year that will be included in
# climatological mean.  (This must be an interval within the 
# observed period).  

# There is a note in the bcgen README saying that the climatology period 
# should generally be left as 1982 - 2001. I'm not sure what 
# would be best to do here, so I will leave it as this.
export mon1clm=1
export iyr1clm=1982
export monnclm=12
export iyrnclm=2001

# export mon1clm=1
# export iyr1clm=1870
# export monnclm=3
# export iyrnclm=2022

# SPECIFY the first and last month and year written to the output 
# file.  (Try to begin at east a few months after the mon1rd,
# iyr1rd, and end a few months before monnrd, iyrnrd, to avoid
# sensitivity to the artificial data outside the observed 
# period.)  

# SW: In bcgen.f90 it's recommended to exclude 2-3 months on either side
# I will exclude first year, as TOGA script requires January start month
export mon1out=1
export iyr1out=1871
export monnout=9
export iyrnout=2016


echo "Initial Checks"
timeinfo=$(ncdump -h  "${dir_data}/$input_file"| grep "float time")
if [ -z "$timeinfo" ]; 
then 
    echo 'Error: Please ensure time variable is named "time" and has type float'
    exit 1
fi 


#######################################################
# copy the target files into a working directory
# Modify these lines for your requirements
#######################################################
echo "organising and renaming data"

cd ${dir_data}
# I'm sort of assuming monnrd is december below with the end day: 
# manually adjust last day of month if necessary

ncks -d time,"${iyr1rd}-${mon1rd}-01 0:00:0.0","${iyrnrd}-${monnrd}-31 00:00:0.0" ${dir_data}/${input_file}  ${dir_temp}/start_file.nc

cd ${dir_temp}

#######################################################
# Make SST file
#######################################################
# SW: they rename sst to SST_cpl
cp start_file.nc temp.nc
# SW: My sst variable is SST
ncrename -v SST,SST_cpl temp.nc sst_cpl.nc
rm temp.nc

#######################################################
# Make sea ice file
#######################################################
cp start_file.nc temp.nc

# SW: here they rename the ice variable
# SW: my ice value is SEAICE
ncrename -v SEAICE,ice_cov temp.nc temp2.nc

# Convert sea ice percentage to fraction
ncap2 -s 'ice_cov=ice_cov/100.' temp2.nc ice_cov.nc
rm temp.nc temp2.nc

#######################################################
# Convert SST file into the standard format
#######################################################
# Modify fill values in the sst_cpl file (which are over land points) to
# have value -1.8 and remove fill and missing value designators;
# change coordinate lengths and names:

# to accomplish this, first run ncdump,
ncdump -d 9,17 sst_cpl.nc > sst_cpl

# then replace _ with -1.8 in SST_cpl,
# then remove lines with _FillValue and missing_value.
sed -i "s/\b_\b/-1.8/g"  sst_cpl
grep -v "_FillValue" sst_cpl > sst_cpl_inter
grep -v "missing_value" sst_cpl_inter > sst_cpl_sub

# To change coordinate lengths and names,
# replace nlon by lon, nlat by lat, TLONG by lon, TLAT by lat.
sed -i "s/\bTLONG\b/lon/g"   sst_cpl_sub
sed -i "s/\bTLON\b/lon/g"   sst_cpl_sub
sed -i "s/\bTLAT\b/lat/g"   sst_cpl_sub
sed -i "s/\bni\b/lon/g"   sst_cpl_sub
sed -i "s/\bnj\b/lat/g"   sst_cpl_sub

# ncdumpThe last step is to run ncgen.
ncgen -o sst_cpl_new.nc sst_cpl_sub
rm sst_cpl sst_cpl_inter sst_cpl_sub

#######################################################
# Convert Sea Ice file into the standard format
#######################################################
# Modify fill values in the ice_cov file (which are over land points)
# to have value 1 and remove fill and missing value designators;
# change coordinate lengths and names; patch longitude and latitude to
# replace missing values: to accomplish this,

# first run ncdump,
ncdump -d 9,17 ice_cov.nc > ice_cov

# then replace _ with 1 in ice_cov,
# then remove lines with _FillValue and missing_value.
sed -i "s/\b_\b/1/g"  ice_cov
grep -v "_FillValue" ice_cov > ice_cov_inter
grep -v "missing_value" ice_cov_inter > ice_cov_sub

# To change coordinate lengths and names,
# replace ni by lon, nj by lat, TLON by lon, TLAT by lat.
sed -i "s/\bTLONG\b/lon/g"  ice_cov_sub
sed -i "s/\bTLON\b/lon/g"   ice_cov_sub
sed -i "s/\bTLAT\b/lat/g"   ice_cov_sub
sed -i "s/\bni\b/lon/g"     ice_cov_sub
sed -i "s/\bnj\b/lat/g"     ice_cov_sub

# To patch longitude and latitude arrays,
# replace values of those arrays with those in sst_cpl file.
# (Note: the replacement of longitude and latitude missing values
# by actual values should not be necessary but is safer.)
# by CD: This step is omitted becuase in the slab ocean run,
# both sst and sea ice are output by the cice module,
# and should have the same longitude and latitude

ncgen -o ice_cov_new.nc ice_cov_sub
rm ice_cov ice_cov_inter ice_cov_sub

#######################################################
# Combining the two files,
# by appending the cice into the sst file
# And Rename the file to ssticetemp.nc.
#######################################################
cp ice_cov_new.nc ice_cov_new2.nc
cp sst_cpl_new.nc sst_cpl_new2.nc
ncks -A -v ice_cov ice_cov_new2.nc sst_cpl_new2.nc
cp sst_cpl_new2.nc ssticetemp.nc

#######################################################
# At this point, you have an SST/ICE file in the correct format.
# However, due to CAM's linear interpolation between mid-month values,
# you need to apply a procedure to assure that the computed monthly means are consistent with the input data.
# To do this, you can invoke the bcgen code in models/atm/cam/tools/icesst and following the following steps:
# reference: https://bb.cgd.ucar.edu/problem-recipe-using-b-compset-output-create-sstice-forcing-files-f-compset
# This toolbox does not come with the CESM but can be found here: ftp://ftp.cgd.ucar.edu/archive/SSTICE/
# Credit goes to Jim Rosinski, Brian Eaton, B. Eaton for coding up this toolbox
# It now comes with the code in Github, with all following modifications:
# 1. In driver.f90, sufficiently expand the lengths of variables prev_history and history
#    (16384 should be sufficient); also comment out the test that the
#    climate year be between 1982 and 2001 (lines 152-158).
# 2. In bcgen.f90 and setup_outfile.f90, change the dimensions of xlon and xlat to (nlon,nlat);
#    this is to accommodate use of non-cartesian ocean grid.
# 3. In setup_outfile.f90, modify the 4th and 5th arguments in the calls to
#    wrap_nf_def_var for lon and lat to be 2 and dimids;
#    this is to accommodate use of non-cartesian ocean grid.
# 4. Adjust Makefile to have proper path for LIB_NETCDF and INC_NETCDF.
#######################################################

# Rename SST_cpl to SST, and ice_cov to ICEFRAC in the current SST/ICE file:
cp ssticetemp.nc ssticetemp_new.nc
ncrename -v SST_cpl,SST -v ice_cov,ICEFRAC ssticetemp_new.nc

# Modify namelist accordingly.
cd ${dir_tool}
cd bcgen

sed -i "s/\( mon1rd = \).*/\1${mon1rd}/"                 namelist
sed -i "s/\( iyr1rd = \).*/\1${iyr1rd}/"                 namelist
sed -i "s/\( monnrd = \).*/\1${monnrd}/"                 namelist
sed -i "s/\( iyrnrd = \).*/\1${iyrnrd}/"                 namelist

sed -i "s/\( mon1 = \).*/\1${mon1}/"                     namelist
sed -i "s/\( iyr1 = \).*/\1${iyr1}/"                     namelist
sed -i "s/\( monn = \).*/\1${monn}/"                     namelist
sed -i "s/\( iyrn = \).*/\1${iyrn}/"                     namelist

sed -i "s/\( mon1clm = \).*/\1${mon1clm}/"               namelist
sed -i "s/\( iyr1clm = \).*/\1${iyr1clm}/"               namelist
sed -i "s/\( monnclm = \).*/\1${monnclm}/"               namelist
sed -i "s/\( iyrnclm = \).*/\1${iyrnclm}/"               namelist

sed -i "s/\( mon1out = \).*/\1${mon1out}/"               namelist
sed -i "s/\( iyr1out = \).*/\1${iyr1out}/"                 namelist
sed -i "s/\( monnout = \).*/\1${monnout}/"               namelist
sed -i "s/\( iyrnout = \).*/\1${iyrnout}/"                 namelist



echo "compiling"

# Make bcgen and execute per instructions.
gmake

# Run the bcgen code and the resulting sstice_ts.nc
# file is the desired ICE/SST file.
rm  ssticetemp_new.nc
rm  *.nc
ln -sf ${dir_temp}/ssticetemp_new.nc .

echo "running bcgen"

./bcgen -i ssticetemp_new.nc -c sstice_clim.nc -t sstice_timeseries.nc < namelist

# Place new SST/ICE file in desired location.
cp -i sstice_clim.nc ${dir_output}
cp -i sstice_timeseries.nc ${dir_output}

# SW: remove the outputs from the bcgen directory
rm sstice_clim.nc
rm sstice_timeseries.nc

# SW: empty the Temp directory
rm ${dir_temp}/*