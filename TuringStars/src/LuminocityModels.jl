module LuminocityModels

using Turing
using StructTypes
using Meshes
using Revise

include("Roche.jl")
include("LuminocityFunctions.jl")

using .Roche
using .LuminocityFunctions

export
    MeshParams,
    ChannelParams,
    ModelParams,
    ChainParams,
    zeroth_model,
    first_model_legacy,
    mq_uniform_model,
    mq_inverted_model,
    star_magnitude,
    star_magnitude_smoker,
    @dist_or_const


@kwdef struct MeshParams
    catmullclark_iterations::Int = 4
    mass_quotient_nodes = 0.01 : 0.01 : 10.
end

@kwdef struct ModelParams
    model_function::Function = zeroth_model
    period::Union{Float64, Distribution}
    β::Union{Float64, Distribution} = 0.08
    temperature_at_bottom::Union{Float64, Distribution} = 3500.
    m_dwarf::Distribution = Uniform(0.3, 1.44)
    m_giant::Distribution = Uniform(0.6, 10.0)
    mass_quotient::Distribution = Uniform(0.1, 10.0)
    mass_quotient_inv::Distribution = Uniform(0.1, 10.0)
end

StructTypes.StructType(::Distribution) = StructTypes.StringType()

@kwdef struct ChannelParams
    measurements_t::Vector{Float64}
    measurements_y::Vector{Float64}
    darkening_function = one1
    darkening_coefs_interpolant = T -> ()
    luminocity_function = T_4
    σ_measured::Vector{Float64}
end
@kwdef struct ChainParams
    mesh_params::MeshParams = MeshParams()
    model_params::ModelParams = ModelParams()
    channels::Vector{ChannelParams}
    n_samples::Int
    n_chains::Int = 1
    sampler::Turing.InferenceAlgorithm = NUTS()
    init_params::Union{Nothing, NamedTuple} = nothing
end
StructTypes.StructType(::Turing.InferenceAlgorithm) = StructTypes.StringType()



macro dist_or_const(var, val)
    esc(:(
        if isa($val, Distribution)
            $var ~ $val
        elseif isa($val, Float64) || isa($val, Vector{Float64})
            $var = $val
        else
            throw(ArgumentError("val must be Distribution, Float64 or Vector{Float64}"))
        end
    ))
end


function star_magnitude(phases; mass_quotient, observer_angle,
                        temperature_at_bottom, β, interpolated_mesh,
                        luminocity_function, darkening_function, darkening_coefs_interpolant)

    directions = [(
        sin(observer_angle) * cos(phase),
        sin(observer_angle) * sin(phase),
        cos(observer_angle)
    ) for phase ∈ phases]

    function temperature(g)
        return temperature_at_bottom * abs(g)^β
    end

    mesh = interpolated_mesh(mass_quotient)
    mesh = avg_over_faces(mesh, :g)
    mesh = apply_function(mesh, temperature, :g, :T)
    mesh = apply_function(mesh, luminocity_function, :T, :L)
    mesh = apply_function(mesh, darkening_coefs_interpolant, :T, :darkening_coefs)

    normals = calc_function_on_faces(mesh, normalized_normal)
    areas = calc_function_on_faces(mesh, area)

    luminocities = [
        integrate_data_over_mesh(mesh, :L, direction, normals, areas, darkening_function)
        for direction ∈ directions
    ]
    return @. -2.5 * log10(luminocities)
end

star_magnitude(
    phases,
    interpolated_mesh,
    model_params::ModelParams,
    channel::ChannelParams,
    sample
) = star_magnitude(
    phases;
    mass_quotient = sample[:mass_quotient],
    observer_angle = sample[:observer_angle],
    model_params.temperature_at_bottom,
    model_params.β,
    interpolated_mesh,
    channel.luminocity_function,
    channel.darkening_function,
    channel.darkening_coefs_interpolant
)

star_magnitude_smoker(
    phases,
    interpolated_mesh,
    model_params::ModelParams,
    channel::ChannelParams,
    sample
) = star_magnitude(
    phases;
    mass_quotient = 1 ./ sample[:mass_quotient],
    observer_angle = sample[:observer_angle],
    model_params.temperature_at_bottom,
    model_params.β,
    interpolated_mesh,
    channel.luminocity_function,
    channel.darkening_function,
    channel.darkening_coefs_interpolant
)


function zeroth_model(mesh_params, model_params, channels)
    interpolated_mesh = InterpolatedRocheMesh(mesh_params)

    @model function model(channels::Vector{ChannelParams}, measurements_y = Float64[])
        m_giant ~ model_params.m_giant
        m_dwarf ~ model_params.m_dwarf
        mass_quotient := m_dwarf / m_giant
        mass_quotient_inv := 1. / mass_quotient

        cos_i ~ Uniform(0., 1.)
        observer_angle := acos(cos_i)

        @dist_or_const temperature_at_bottom model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        offset ~ filldist(Flat(), length(channels))
        log_σ_common ~ filldist(Flat(), length(channels))
        σ_common = exp.(log_σ_common)

        for (i, channel) ∈ enumerate(channels)
            phases = initial_phase .+ channel.measurements_t * 2π / period

            predicted_magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, (; mass_quotient, observer_angle))

            predicted_magnitudes .+= offset[i]

            σ = @. √(channel.σ_measured^2 + σ_common[i]^2)

            measurements_y = channel.measurements_y
            measurements_y .~ Normal.(predicted_magnitudes, σ)
        end
    end

    return model(channels)
end

StructTypes.StructType(::typeof(zeroth_model)) = StructTypes.StringType()


function first_model_legacy(mesh_params, model_params, channels)
    interpolated_mesh = InterpolatedRocheMesh(mesh_params)

    @model function model(channels::Vector{ChannelParams}, measurements_y = Float64[])
        mass_quotient ~ model_params.mass_quotient
        mass_quotient_inv := 1. / mass_quotient
        observer_angle ~ Uniform(0., π/2)

        @dist_or_const temperature_at_bottom model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        offset ~ filldist(Flat(), length(channels))
        σ_common = filldist(FlatPos(0.), length(channels))

        for (i, channel) ∈ enumerate(channels)
            phases = initial_phase .+ channel.measurements_t * 2π / period

            predicted_magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, (; mass_quotient, observer_angle))

            predicted_magnitudes .+= offset[i]

            σ = @. √(channel.σ_measured^2 + σ_common[i]^2)

            measurements_y = channel.measurements_y
            measurements_y .~ Normal.(predicted_magnitudes, σ)
        end
    end

    return model(channels)
end

StructTypes.StructType(::typeof(first_model_legacy)) = StructTypes.StringType()



function mq_uniform_model(mesh_params, model_params, channels)
    interpolated_mesh = InterpolatedRocheMesh(mesh_params)

    @model function model(channels::Vector{ChannelParams}, measurements_y = Float64[])
        mass_quotient ~ model_params.mass_quotient
        mass_quotient_inv := 1. / mass_quotient

        cos_i ~ Uniform(0., 1.)
        observer_angle := acos(cos_i)

        @dist_or_const temperature_at_bottom model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        offset ~ filldist(Flat(), length(channels))
        log_σ_common ~ filldist(Flat(), length(channels))
        σ_common = exp.(log_σ_common)

        for (i, channel) ∈ enumerate(channels)
            phases = initial_phase .+ channel.measurements_t * 2π / period

            predicted_magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, (;mass_quotient, observer_angle))

            predicted_magnitudes .+= offset[i]

            σ = @. √(channel.σ_measured^2 + σ_common[i]^2)

            measurements_y = channel.measurements_y
            measurements_y .~ Normal.(predicted_magnitudes, σ)
        end
    end

    return model(channels)
end

StructTypes.StructType(::typeof(mq_uniform_model)) = StructTypes.StringType()



function mq_inverted_model(mesh_params, model_params, channels)
    interpolated_mesh = InterpolatedRocheMesh(mesh_params)

    @model function model(channels::Vector{ChannelParams}, measurements_y = Float64[])
        mass_quotient_inv ~ model_params.mass_quotient_inv
        mass_quotient := 1. / mass_quotient_inv

        cos_i ~ Uniform(0., 1.)
        observer_angle := acos(cos_i)

        @dist_or_const temperature_at_bottom model_params.temperature_at_bottom
        @dist_or_const β model_params.β

        initial_phase ~ Uniform(-π, π)

        @dist_or_const period model_params.period

        offset ~ filldist(Flat(), length(channels))
        log_σ_common ~ filldist(Flat(), length(channels))
        σ_common = exp.(log_σ_common)

        for (i, channel) ∈ enumerate(channels)
            phases = initial_phase .+ channel.measurements_t * 2π / period

            predicted_magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, (;mass_quotient, observer_angle))

            predicted_magnitudes .+= offset[i]

            σ = @. √(channel.σ_measured^2 + σ_common[i]^2)

            measurements_y = channel.measurements_y
            measurements_y .~ Normal.(predicted_magnitudes, σ)
        end
    end

    return model(channels)
end

StructTypes.StructType(::typeof(mq_inverted_model)) = StructTypes.StringType()



end