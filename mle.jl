### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 11decfb6-8f0a-11ee-1bf8-d3faf8756b8b
begin
	import Pkg
	Pkg.activate(".")

	using Revise
	using TuringStars

	using DelimitedFiles
	using DataFrames

	using LombScargle

	using Turing

	using Plots
	using StatsPlots

	using Optim
end

# ╔═╡ 2408848a-e855-44cb-9b4f-3c6ced602e07
using Dates

# ╔═╡ 7feb762a-96d3-4b74-8a0a-cd3ff0c2a352
using HypothesisTests

# ╔═╡ 512ace6e-9008-4d4a-b513-3be764fe3ad4
using LaTeXStrings

# ╔═╡ 6b6cf7e6-340f-4e34-aeef-c7c9f1695380
begin
	plotlyjs()
	theme(:juno)
end

# ╔═╡ 0a209f5e-bd22-4072-9247-c0ddcd42fdb8
begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

# ╔═╡ 67b9b67f-7b71-4e16-904f-c89e244bea7e
period = 227.5687

# ╔═╡ a6dccb45-fa0c-4aff-ae2e-23efc08fd8da
interpolated_mesh = InterpolatedRocheMesh(tetra_sphere(4), 0.1:0.1:10)

# ╔═╡ e36d9949-98cc-4635-835c-400a584fcde4
function group_symbols(sample)
	if isa(sample, Chains)
		sample = get_params(sample)
	end
	sample
end

# ╔═╡ a10f8c4b-7d02-43d4-b5cf-e869e794c3d6
function plot_template(xlabel = "Julian day % period"; kwargs...)
	plot(
		layout = (2, 1),
		title = ["Спектральный канал K" "Спектральный канал J"],
		legend = false,
		xlabel = ["" xlabel],
		ylabel = "Звездная величина",
		yflip = true,
		size = (600, 600),
		margin = 12Plots.px,
	)
	plot!(; kwargs...)
end

# ╔═╡ 471730f6-cd54-45cb-871d-dfdc7d7e4d65
function plot_points_days(channels, period, p = nothing; kwargs...)
	if p == nothing
		p = plot_template()
	end

	for (channel, subplot) ∈ zip(channels, p.subplots)
		scatter!(
			subplot,
			channel.measurements_t .% period,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2,
		)
	end

	scatter!(; kwargs...)
end

# ╔═╡ 16711745-a609-4add-8118-3918af8293aa
function plot_lines_days(model_params, interpolated_mesh,
								channels, samples, p = nothing; kwargs...)
	if p == nothing
		p = plot_template()
	end

	days = 0 : period
	phases_ = @. (days % model_params.period) / model_params.period * 2π

	for (c, (channel, subplot)) ∈ enumerate(zip(channels, p.subplots))
		for i ∈ 1 : length(samples)
			sample = group_symbols(samples[i])

			phases = reshape(phases_ .+ sample[:initial_phase], :)
			magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]

			plot!(
				subplot,
				days .% model_params.period,
				magnitudes;
				kwargs...
			)
		end
	end
	p
end

# ╔═╡ 82133dc0-bdd2-4df9-83df-3e9ccb9d3df1
function plot_points_phases(channels, period, initial_phase, p = nothing; kwargs...)
	if p == nothing
		p = plot_template("Фаза")
	end

	for (channel, subplot) ∈ zip(channels, p.subplots)
		phases = (channel.measurements_t .% period) ./ period .+ (initial_phase / 2π)
	
		scatter!(
			subplot,
			phases,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2;
			kwargs...
		)
		scatter!(
			subplot,
			phases .- 1,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2;
			kwargs...
		)
		scatter!(
			subplot,
			phases .+ 1,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2;
			kwargs...
		)
	end
	scatter!(xlim = (-0.25, 1.25))
end

# ╔═╡ 3460e194-76f9-4c06-a7a1-b85845b8ea9f
function plot_lines_phases(model_params, interpolated_mesh,
							channels, samples, p = nothing; kwargs...)
	if p == nothing
		p = plot_template("Фаза")
	end

	phases = -0.25 : 0.01 : 1.25

	for (c, (channel, subplot)) ∈ enumerate(zip(channels, p.subplots))

		for s ∈ 1 : length(samples)
			sample = group_symbols(samples[s])

			vals = star_magnitude(phases .* 2π, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]

		plot!(
			subplot,
			phases,
			vals;
			kwargs...
		)	
		end
	end

	return p
end

# ╔═╡ 7773cf80-05cc-4c4d-9185-a81b1a237fa1
model_0 = ModelParams(; period);

# ╔═╡ 15fa2682-e8e7-44b2-993a-3d1c7f66cf4e
model_mq_uniform = ModelParams(; period, model_function = mq_uniform_model);

# ╔═╡ fca42ad4-ef24-4480-8573-6f5b3dab8584
model_mq_inverted = ModelParams(; period, model_function = mq_inverted_model);

# ╔═╡ f3e28f2a-c23e-42d1-9b8a-467270725aa8
model_bb = ModelParams(; period, model_function = mq_inverted_model);

# ╔═╡ 5449571a-5e90-4686-84de-6bfe09fba48f
channels_0 = [
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.K,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = K_coefs_interpolant,
		luminocity_function = phoenixK,
		σ_measured = points.K_err,
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = J_coefs_interpolant,
		luminocity_function = phoenixJ,
		σ_measured = points.J_err,
	)
];

# ╔═╡ d2bf9a47-5c7a-4db8-b028-8d3f7383ff95
instance_0 = model_0.model_function(MeshParams(), model_0, channels_0);

# ╔═╡ 4dd7d6f5-0711-4156-b2e1-20a017d5276c
instance_mq_uniform = model_mq_uniform.model_function(MeshParams(), model_mq_uniform, channels_0);

# ╔═╡ 351a6436-557d-4d76-871e-1f9e466de676
instance_mq_inv = model_mq_inverted.model_function(MeshParams(), model_mq_inverted, channels_0);

# ╔═╡ daf42ed3-feaa-45ea-969e-ac91f1f3ea5b
channels_bb = [
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.K,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = K_coefs_interpolant,
		luminocity_function = black_body_K,
		σ_measured = points.K_err,
	)
	ChannelParams(
		measurements_t = points.day,
		measurements_y = points.J,
		darkening_function = claret_darkening,
		darkening_coefs_interpolant = J_coefs_interpolant,
		luminocity_function = black_body_J,
		σ_measured = points.J_err,
	)
];

# ╔═╡ 49f0671b-3110-48bd-a8f3-5813a5634af7
instance_bb = model_bb.model_function(MeshParams(), model_bb, channels_bb);

# ╔═╡ f278a1e6-f459-4582-9400-9bfa2d852d14
begin
	local m_giant = 0.6
	local m_dwarf = 1.0
	local mass_quotient = m_dwarf / m_giant
	local mass_quotient_inv = m_giant / m_dwarf

	local cos_i = 0.6
	local observer_angle = acos(cos_i)

	local initial_phase = -1.46867
	local offset = [46.6894, 48.8942]
	local log_σ_common = log.([0.01, 0.01])

	init_params_0 = (; m_giant, m_dwarf, mass_quotient, mass_quotient_inv, cos_i, observer_angle, initial_phase, offset, log_σ_common)
end

# ╔═╡ 535fb7b7-781c-4a48-b5f4-64aa81ff8bed
begin
	local p = plot_points_days(channels_0, period)
	plot_lines_days(model_0, interpolated_mesh, channels_0, [init_params_0], p)
end

# ╔═╡ d7bf9cbb-768e-4425-97be-b2914ef7d72d
begin
	local m_giant = 0.8
	local m_dwarf = 1.0
	local mass_quotient = m_dwarf / m_giant
	local mass_quotient_inv = m_giant / m_dwarf

	local cos_i = 0.5
	local observer_angle = acos(cos_i)

	local initial_phase = -1.46867
	local offset = [18.85, 21.17]
	local log_σ_common = log.([0.01, 0.01])

	init_params_bb = (; m_giant, m_dwarf, mass_quotient, mass_quotient_inv, cos_i, observer_angle, initial_phase, offset, log_σ_common)
end

# ╔═╡ b1701d19-4cb7-4c35-943e-6a9cd415c4a1
function logpdf_objective(model)
	var_info = DynamicPPL.VarInfo(model)

	return function f(vars)
		var_info_ = DynamicPPL.unflatten(var_info, vars)
		return logjoint(model, var_info_)
	end
end

# ╔═╡ a3b38672-7e73-4701-abc6-30d057f09754
function MAP_point(model, initial_params)
	var_info = DynamicPPL.VarInfo(model)
	syms = DynamicPPL.syms(var_info)
	x0 = vcat(collect(initial_params[syms])...)

	optim_result = maximize(logpdf_objective(model), x0, NelderMead())
	optim_array = Optim.maximizer(optim_result)

	var_info = DynamicPPL.unflatten(var_info, optim_array)
	tuples = Turing.Inference.getparams(model, var_info)
	return Chains(reshape(last.(tuples), (1, :, 1)), Symbol.(first.(tuples)))
end

# ╔═╡ f7c074e7-d5b1-4b2c-9ab6-57e85b237a04
MAP_0 = MAP_point(instance_0, init_params_0)

# ╔═╡ 67391953-6f49-466e-bc72-349d9df43b0b
MAP_mq_uniform = MAP_point(instance_mq_uniform, init_params_0)

# ╔═╡ 09c08a94-a567-4ba8-af64-424e9dbbe75b
MAP_mq_inv = MAP_point(instance_mq_inv, init_params_0)

# ╔═╡ c5bc01de-c96c-40d2-9b9f-d068ec054798
MAP_bb = MAP_point(instance_bb, init_params_bb)

# ╔═╡ 8cf3d1d8-503f-47ec-abf1-b8e2293b96bf
# begin
# 	var_info = DynamicPPL.VarInfo(instance_0)
# 	syms = DynamicPPL.syms(var_info)
# 	x0 = vcat(collect(init_params_0[syms])...)
# 	x1 = MAP_0.value[:, [:m_giant, :m_dwarf, :cos_i, :initial_phase, Symbol("offset[1]"), Symbol("offset[2]"), Symbol("log_σ_common[1]"), Symbol("log_σ_common[2]")]]
# 	x1 = reshape(Array(x1), :)
# end

# ╔═╡ 2e887134-51f0-421e-9a89-c939c714a671
# typeof(instance_0)

# ╔═╡ 14c44a4c-d5e0-4e04-a74d-506695b0a499
# MAP_0.value[:, [:m_giant, :m_dwarf, :cos_i, :initial_phase, Symbol("offset[1]"), Symbol("offset[2]"), Symbol("log_σ_common[1]"), Symbol("log_σ_common[2]")]]

# ╔═╡ 19fc43ce-0061-46a3-a5fc-018d1f46a044
# logpdf_objective(instance_0)(x0)

# ╔═╡ 0bac0f09-c2ab-41fb-8866-5695693e4960
# logpdf_objective(instance_0)(x1)

# ╔═╡ d1e094f8-a0fd-45ce-9255-4b9de702433b
begin
	plotlyjs()
	p = plot_points_days(channels_0, period)
	plot_lines_days(model_mq_inverted, interpolated_mesh, channels_0, [MAP_mq_inv], p; label = "q ~ Uniform")
	plot_lines_days(model_mq_uniform, interpolated_mesh, channels_0, [MAP_mq_uniform], p; label = "q^-1 ~ Uniform")
	plot_lines_days(model_0, interpolated_mesh, channels_0, [MAP_0], p; label = "m_giant, m_dwarf ~ Uniform")
	plot_lines_days(model_bb, interpolated_mesh, channels_bb, [MAP_bb], p; label = "blackbody; q ~ Uniform")
	plot!(legend = true)
end

# ╔═╡ 1edd82ba-6906-4ff5-8982-6020dd2b830f
map_plot = begin
	gr()
	theme(:default)
	local p = plot_points_phases(channels_0, period, mean(mean([MAP_mq_inv[:initial_phase], MAP_bb[:initial_phase]])), markercolor = "black", label = "")

	plot_lines_phases(model_mq_inverted, interpolated_mesh, channels_0, [MAP_mq_inv], p; label = "PHOENIX", linecolor = 1)

	plot_lines_phases(model_bb, interpolated_mesh, channels_bb, [MAP_bb], p; label = "black body", linecolor = 2)
 
	plot!(plot_title = "MLE-оценки кривых блеска", plot_titlevspan = 0.065, bottom_margin = 0Plots.px)
	plot!(p.subplots[2], legend = :bottom, legend_title = "Модель излучения", legend_title_font_pointsize = 8)
end

# ╔═╡ b9a45bc0-21b6-4eed-b06f-0893c20480c4
savefig(map_plot, "tex/pic_drafts/mle_different_L.pdf")

# ╔═╡ 66651fe4-7115-40e9-935d-95795757da4e
md"### ``\chi^2``"

# ╔═╡ 28df81cf-e27f-40ee-9ce6-5ff202994dbc
function chi2_value(model_params, channels, MAP_params)
	MAP_params = group_symbols(MAP_params)

	sum(enumerate(channels)) do (i, channel)

		phases = channel.measurements_t ./ model_params.period .* 2π .+ MAP_params[:initial_phase]

		predicted_magnitudes = MAP_params[:offset][i]
		predicted_magnitudes = star_magnitude(
			phases;
			mass_quotient = MAP_params[:mass_quotient],
			observer_angle = MAP_params[:observer_angle],
			temperature_at_bottom = model_params.temperature_at_bottom,
			interpolated_mesh,
			β = model_params.β,
			luminocity_function = channel.luminocity_function,
			darkening_function = channel.darkening_function,
			darkening_coefs_interpolant = channel.darkening_coefs_interpolant
		)
		predicted_magnitudes .+= MAP_params[:offset][i]

		σ² = @. channel.σ_measured^2 + exp(2MAP_params[:log_σ_common][i])

		sum(@. (channel.measurements_y - predicted_magnitudes)^2 / σ²)
	end
end

# ╔═╡ be6c4b94-861e-4907-8a53-dbc3a674664e
chi2_value_ = chi2_value(model_0, channels_0, MAP_0)

# ╔═╡ bce9145d-23f2-4e22-90bc-0c0089502286
Chi2Dist = Chisq(2 * length(points.day) - 5)

# ╔═╡ 363ea8d3-a174-4186-a967-be659aaa9a72
1 - cdf(Chi2Dist, chi2_value_)

# ╔═╡ eefb61ae-924d-4bb7-9b9a-7235b1f4b826
md"### err = ``f(t)``"

# ╔═╡ b85b3eb7-5257-43d9-9261-0736f7102474
function differences(model_params, channels, MAP_params)
	MAP_params = group_symbols(MAP_params)

	map(enumerate(channels)) do (i, channel)

		phases = channel.measurements_t ./ model_params.period .* 2π .+ MAP_params[:initial_phase]
		phases = reshape(phases, :)

		predicted_magnitudes = MAP_params[:offset][i]
		predicted_magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, MAP_params)
		predicted_magnitudes .+= MAP_params[:offset][i]

		@. channel.measurements_y - predicted_magnitudes
	end
end

# ╔═╡ c5ad427c-c144-4cc2-b531-ff830beaa21e
scatter(
	# julian2datetime.(2450199.5 .+ points.day),
	#(points.day .% estimated_period) / estimated_period,
	sign.(differences(model_0, channels_0, MAP_0)[1]),
	#label = ["K" "J"],
	legend = false,
	markersize = 2
)

# ╔═╡ 6aded0f0-c6b3-4bf2-872c-eb9b264caf00
histogram(differences(model_0, channels_0, MAP_0)[1], nbins = 20)

# ╔═╡ fac51a59-200f-4413-a63e-46571e31f324
md"### Runs Test"

# ╔═╡ e104d63f-603c-452d-916b-08be782834ce
WaldWolfowitzTest(differences(model_0, channels_0, MAP_0)[1])

# ╔═╡ 47d295b7-fb33-40e3-84e2-7d788b3c3acf
function plot_positive_negative(x, y; p = nothing, kwargs...)
	positive_mask = y .> 0
	negative_mask = .! positive_mask

	if p == nothing
		p = plot()
	end

	scatter!(x[negative_mask], y[negative_mask], markercolor = 1; kwargs...)
	scatter!(x[positive_mask], y[positive_mask], markercolor = 2; kwargs...)
end

# ╔═╡ 9bbfe3ec-9f74-4eab-8ea5-57b434abdcaa
plot_positive_negative(y; p = nothing, kwargs...) = plot_positive_negative(1 : length(y), y; p, kwargs...)

# ╔═╡ de5cfa4f-3c91-4522-b01a-3fe4697cff72
gr()

# ╔═╡ ed305c32-10ce-4aa0-9df4-d2bb644c2341
diffs_plot = begin
	gr()
	diffs_JK = differences(model_mq_uniform, channels_0, MAP_mq_uniform)
	diffs_plot = plot_positive_negative(diffs_JK[1]; markershape = :+, label = false)
	plot_positive_negative(diffs_JK[2]; p = diffs_plot, markersize = 2, markerstrokewidth = 0, label = false)
	plot!(
		size = (600, 350),
		xlabel = "Номер точки",
		ylabel = L"m_\mathrm{meas} - m_\mathrm{model}",
		title = "Разность между измерениями\nи модельной кривой",
		top_margin = 30Plots.px
	)
	scatter!([], markershape = :+, markersize = 3, markercolor = "grey", label = "Канал K")
	scatter!([], markersize = 1.5, markercolor = "grey", markerstrokecolor = "grey", label = "Канал J")
end

# ╔═╡ 992778bc-59bd-40e1-b313-66c1f2ef34c3
savefig(diffs_plot, "tex/pic_drafts/residuals.pdf")

# ╔═╡ Cell order:
# ╠═11decfb6-8f0a-11ee-1bf8-d3faf8756b8b
# ╠═6b6cf7e6-340f-4e34-aeef-c7c9f1695380
# ╠═0a209f5e-bd22-4072-9247-c0ddcd42fdb8
# ╠═67b9b67f-7b71-4e16-904f-c89e244bea7e
# ╠═a6dccb45-fa0c-4aff-ae2e-23efc08fd8da
# ╠═e36d9949-98cc-4635-835c-400a584fcde4
# ╠═a10f8c4b-7d02-43d4-b5cf-e869e794c3d6
# ╠═471730f6-cd54-45cb-871d-dfdc7d7e4d65
# ╠═16711745-a609-4add-8118-3918af8293aa
# ╠═82133dc0-bdd2-4df9-83df-3e9ccb9d3df1
# ╠═3460e194-76f9-4c06-a7a1-b85845b8ea9f
# ╠═7773cf80-05cc-4c4d-9185-a81b1a237fa1
# ╠═15fa2682-e8e7-44b2-993a-3d1c7f66cf4e
# ╠═fca42ad4-ef24-4480-8573-6f5b3dab8584
# ╠═f3e28f2a-c23e-42d1-9b8a-467270725aa8
# ╠═d2bf9a47-5c7a-4db8-b028-8d3f7383ff95
# ╠═4dd7d6f5-0711-4156-b2e1-20a017d5276c
# ╠═351a6436-557d-4d76-871e-1f9e466de676
# ╠═49f0671b-3110-48bd-a8f3-5813a5634af7
# ╠═5449571a-5e90-4686-84de-6bfe09fba48f
# ╠═daf42ed3-feaa-45ea-969e-ac91f1f3ea5b
# ╠═f278a1e6-f459-4582-9400-9bfa2d852d14
# ╠═535fb7b7-781c-4a48-b5f4-64aa81ff8bed
# ╠═d7bf9cbb-768e-4425-97be-b2914ef7d72d
# ╠═b1701d19-4cb7-4c35-943e-6a9cd415c4a1
# ╠═a3b38672-7e73-4701-abc6-30d057f09754
# ╠═f7c074e7-d5b1-4b2c-9ab6-57e85b237a04
# ╠═67391953-6f49-466e-bc72-349d9df43b0b
# ╠═09c08a94-a567-4ba8-af64-424e9dbbe75b
# ╠═c5bc01de-c96c-40d2-9b9f-d068ec054798
# ╠═8cf3d1d8-503f-47ec-abf1-b8e2293b96bf
# ╠═2e887134-51f0-421e-9a89-c939c714a671
# ╠═14c44a4c-d5e0-4e04-a74d-506695b0a499
# ╠═19fc43ce-0061-46a3-a5fc-018d1f46a044
# ╠═0bac0f09-c2ab-41fb-8866-5695693e4960
# ╠═d1e094f8-a0fd-45ce-9255-4b9de702433b
# ╠═1edd82ba-6906-4ff5-8982-6020dd2b830f
# ╠═b9a45bc0-21b6-4eed-b06f-0893c20480c4
# ╟─66651fe4-7115-40e9-935d-95795757da4e
# ╠═28df81cf-e27f-40ee-9ce6-5ff202994dbc
# ╠═be6c4b94-861e-4907-8a53-dbc3a674664e
# ╠═bce9145d-23f2-4e22-90bc-0c0089502286
# ╠═363ea8d3-a174-4186-a967-be659aaa9a72
# ╟─eefb61ae-924d-4bb7-9b9a-7235b1f4b826
# ╠═2408848a-e855-44cb-9b4f-3c6ced602e07
# ╠═b85b3eb7-5257-43d9-9261-0736f7102474
# ╠═c5ad427c-c144-4cc2-b531-ff830beaa21e
# ╠═6aded0f0-c6b3-4bf2-872c-eb9b264caf00
# ╟─fac51a59-200f-4413-a63e-46571e31f324
# ╠═7feb762a-96d3-4b74-8a0a-cd3ff0c2a352
# ╠═e104d63f-603c-452d-916b-08be782834ce
# ╠═47d295b7-fb33-40e3-84e2-7d788b3c3acf
# ╠═9bbfe3ec-9f74-4eab-8ea5-57b434abdcaa
# ╠═512ace6e-9008-4d4a-b513-3be764fe3ad4
# ╠═de5cfa4f-3c91-4522-b01a-3fe4697cff72
# ╠═ed305c32-10ce-4aa0-9df4-d2bb644c2341
# ╠═992778bc-59bd-40e1-b313-66c1f2ef34c3
