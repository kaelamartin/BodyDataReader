module BodyDataReader

using JLD2
using FileIO
using FTPClient

using Printf

using Dates
using DelimitedFiles
using LinearAlgebra
using Statistics
using SparseArrays

##########
# Boddat #
##########

include("boddat.jl")

export boddat,
	   boddatephem!,
	   boddatorient

##############
# Dictionary #
##############

include("dict.jl")

export gendict!,
	   savedict,
	   tryloaddict!,
	   MB,
	   SB

###########
# General #
###########

include("general.jl")

export getspk,
	   getname,
	   getgm,
	   getj2,
	   getradius,
	   gettriaxial,
	   getrotperiod,
	   getcb

##############
# Small Body #
##############

include("smallbody.jl")

export getmagnitude,
	   getocc,
	   gettype,
	   getclose

#############
# Ephemeris #
#############

include("ephemeris.jl")

export getreference,
	   getdate,
	   ephem!,
	   ephem,
	   EphemType,
	   TEPHF

###############
# Orientation #
###############

include("orientation.jl")

export pole,
	   poledcm,
	   pm,
	   pmdcm,
	   eqx,
	   eqxdcm,
	   qpole,
	   qpoledcm,
	   qpm,
	   qpmdcm,
	   qeqx,
	   qeqxdcm

##########
# Spline #
##########

include("spline.jl")

end
