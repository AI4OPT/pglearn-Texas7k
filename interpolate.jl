using Base.Threads
using Dates
using Printf

using HDF5
using Interpolations

"""
    _check_texas7k_demand_data(D; check_type::Bool=true, check_size::Bool=true)
# Arguments
* `D::Dict`

# Keyword arguments
* `check_type::Bool=true`: whether to check value types
* `check_type::Bool=true`: whether to check value size

# Throws
An error is thrown if some data checks fail
"""
function _check_texas7k_demand_data(D; check_type::Bool=true, check_size::Bool=true)

    TEXAS7K_DEMAND_SCHEMA = Dict(
        "datetime" => Dict(
            "type" => Vector{String},
            "size" => (8760,),
        ),
        "pd" => Dict(
            "type" => Matrix{Float32},
            "size" => (8760, 4549),
        ),
        "qd" => Dict(
            "type" => Matrix{Float32},
            "size" => (8760, 4549),
        ),
    )

    for (k, s) in TEXAS7K_DEMAND_SCHEMA
        haskey(D, k) || error("Missing key: $k")
        # check type
        v = D[k]
        (check_type && isa(v, s["type"])) || error(
            """Invalid value type for key $k:
            * Expected: $(s["type"])
            * Actual  : $(typeof(v))"""
        )
        (check_size && size(v) == s["size"]) || error(
            """Invalid size of data for key $k:
            * Expected: $(s["size"])
            * Actual  : $(size(v))"""
        )
    end

    return nothing
end

function _load_texas7k_demand_2020_consolidated()
    fpath = joinpath(@__DIR__, "data", "texas7k_demand_2020.h5")
    if !isfile(fpath)
        error("Missing file: ", abspath(normpath(fpath)))
    end

    D = h5read(joinpath(@__DIR__, "data", "texas7k_demand_2020.h5"), "/")

    _check_texas7k_demand_data(D)
    return D
end

function _load_texas7k_demand_2020_monthly()
    # Check that all files exist
    fpaths = [
        joinpath(@__DIR__, "data", @sprintf("texas7k_demand_2020-%02d.h5", mm))
        for mm in 1:12
    ]
    all(isfile.(fpaths)) || error("Some monthly h5 files are missing; please check that you cloned the repository correctly")
    
    # All files exist, let's load them and consolidate into a single dictionary
    ds = [h5read(fpath, "/") for fpath in fpaths]
    D = Dict(
        k => reduce(vcat, [d[k] for d in ds])
        for k in ["datetime", "pd", "qd"]
    )

    _check_texas7k_demand_data(D)
    return D
end

"""
    load_texax7k_demand_2020()

Load the Texas7k 2020 demand time series data.

Returns a Dictionary with the following keys
* `datetime::Vector{String}`
* `pd::Matrix{Float32}`: a 8760×4549 Matrix of hourly nodal active power demand
* `qd::Matrix{Float32}`: a 8760×4549 Matrix of hourly nodal reactive power demand
"""
function load_texas7k_demand_2020()
    if isfile(joinpath(@__DIR__, "data", "texas7k_demand_2020.h5"))
        return _load_texas7k_demand_2020_consolidated()
    else
        return _load_texas7k_demand_2020_monthly()
    end
end

"""
    interpolate_to_5min(D)

Interpolate active and reactive demand from hourly to 5-minute granularity.
The interpolation is done using cubic splines.

# Arguments
* `D::Dict`: dictionary containing 3 keys:
    * `datetime`
    * `D["pd"]` is a Matrix{Float32} of size `T*L`, such that 
        D["pd"][t, i] is the active demand of load `i` at time step `t`.
    * `D["qd"]` is a Matrix{Float32} of size `T*L`, such that 
        D["qd"][t, i] is the active demand of load `i` at time step `t`.

# Returns
* A dictionary with same keys as `D`, but with 5-min interpolated data
"""
function interpolate_to_5min(D)
    # Hourly data info
    dts_1hr = DateTime.(D["datetime"])
    T_1hr = length(dts_1hr)
    pd_1hr = D["pd"]
    qd_1hr = D["qd"]
    L = size(pd_1hr, 2)  # number of loads

    dts_5min = collect(minimum(dts_1hr):Minute(5):maximum(dts_1hr))
    # The interpolation will create an additional 11 points per hour,
    #   except for the last hour
    T_5min = 12 * (T_1hr - 1) + 1
    # pre-allocate interpolated active/reactive demand
    pd_5min = zeros(Float32, T_5min, L)
    qd_5min = zeros(Float32, T_5min, L)

    # Now we do the actual interpolation
    for (x_1hr, x_5min) in zip([pd_1hr, qd_1hr], [pd_5min, qd_5min])
        @threads for i in 1:L
            local itp = cubic_spline_interpolation(1:T_1hr, x_1hr[:, i])
            x_5min[:, i] .= itp.(collect(1:(1/12):T_1hr))
        end
    end

    D_5min = Dict(
        "datetime" => string.(dts_5min),
        "pd" => pd_5min,
        "qd" => qd_5min,
    )

    return D_5min
end

function main_interpolate()
    D = load_texas7k_demand_2020()
    D_5min = interpolate_to_5min(D)

    # Save to h5 file
    if !isdir(joinpath(@__DIR__, "data", "5min"))
        @warn "Path `data/5min` does not exist; creating it"
        mkpath(joinpath(@__DIR__, "data", "5min"))
    end

    h5open(joinpath(@__DIR__, "data", "5min", "texas7k_demand_2020_5min.h5"), "w") do fid
        for (k, v) in D_5min
            fid[k] = v
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_interpolate()
    exit(0)
end
