# Table of Contents

- [Introduction](#introduction)
- [Using `boddat`](#using-boddat)
- [Direct Function Calls](#direct-function-calls)
- [Notes](#notes)

# Introduction

BodyDataReader provides two ways to retrieve desired data 

`boddat` is the easiest way to retrieve data using BodyDataReader. Generally, `boddat` is suitable for use in code without great performance concerns and retrieving data of interest to the user from the REPL.

Using `boddat` comes with significant overhead and thus it is recommended to use direct function calls in performance-centric applications. 
However, there are tasks that `boddat` automatically performs that the user must instead perform when calling functions directly.

# Using `boddat`

At minimum, `boddat` takes three arguments:

```julia
boddat("command string", [body], dict)
```

  * `"command string"` specifies which data type to retrieve, listed in the `boddat` help output.
  * `[body]` specifies the body to retrieve data for.
  * `dict` is a dictionary that the program caches data in. It must be of type Dict{String,Any}.
  
When a call is made to `boddat`, it does several of things:

  1. Checks if dict is empty, and tries to load the previous dict from file if it is.
  2. Looks for time input, and uses the current time if none is given.
  3. Parses options if given.
  4. Parses command string.
  5. Calls specified function.
  6. Saves dict with new data to file, if saving is enabled.
  7. Returns the requested data.
  
Every one of these tasks except 5 are superflous to the task of retrieving the requested data. The user would otherwise have to perform these tasks, which is why they are included in `boddat`.

# Direct Function Calls

Direct function calls are far more efficient for retrieving data compared to using `boddat` because they eliminate overhead, or allow the user to batch significant amounts of overhead.

The primary task that the user must perform when using direct function calls compared to `boddat` is managing the dictionary directly. Generally, use of BodyDataReader without `boddat` follows this flow:

  1. A new dictionary is created and populated it either by loading the old dictionary into it or populating it empty.
  2. Various calls to BodyDataReader functions are made.
  3. The dictionary is saved for later sessions.
  
BodyDataReader provides three functions for manipulating dictionaries to use with functions: `gendict!`, `tryloaddict!`, and `savedict`.

`gendict!` populates an empty dictionary with the required structures for use with BodyDataReader functions.

`tryloaddict!` attempts to load the saved dictionary from file, if it exists. It can take an empty dictionary.

`savedict` saves the dictionary to file for later use.

Generally, the use of these commands will follow a typical structure:

```julia
dict = Dict{String,Any}()
tryloaddict!(dict)
...
[function calls]
...
savedict(dict)
```

There are two types of functions for data retrieval: time independent, and time dependent.

Time independent commands generally have the following form:

```julia
function(body, dict)
```

Time dependent commands generally have the following form:

```julia
function(body, time or time range, dict)
```

Some functions, notably ephemeris commands, have additional arguments that are required. Consult the [reference](doc/Reference.md) for information on available functions.

# Notes

Some behaviors may not be entirely obvious or intuitive, so they are explained here.

Data retrieval functions only take SPK IDs as input. The function `getspk` can be used to attempt to parse a string to an SPK ID if unknown.

All times are relative J2000 epoch (Julian date 2451545.0).

All ephemeris data is retrieved relative to the Solar System Barycenter by default. Overloaded functions are available to retrieve ephemeris relative to other bodies.