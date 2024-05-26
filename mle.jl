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

# ╔═╡ 18694e43-32c0-424f-ade0-63c6fc8af923
pgram = lombscargle(points.day, points.K, samples_per_peak = 100)

# ╔═╡ 4b30de9a-e872-422c-9b12-2cc67942aba3
pgramplot = plot(
	freq(pgram),
	power(pgram),
	xlabel = "частота (1/день)",
	title = "Периодограмма Ломба-Скаргла",
	legend = nothing,
	xlims = (0.003, 0.015),
	linewidth = 1.2,
	size = (600, 360)
)

# ╔═╡ 49216477-bda7-487c-9c85-56e1ed842d84
savefig(pgramplot, "tex/pic_drafts/periodogram.svg")

# ╔═╡ c73161c8-0eb3-40af-b527-40a05e10ae24
2findmaxperiod(pgram)[1]

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

# ╔═╡ 7773cf80-05cc-4c4d-9185-a81b1a237fa1
model_0 = ModelParams(; period);

# ╔═╡ 15fa2682-e8e7-44b2-993a-3d1c7f66cf4e
model_q_uniform = ModelParams(; period, model_function = q_uniform_model);

# ╔═╡ fca42ad4-ef24-4480-8573-6f5b3dab8584
model_q_inverted = ModelParams(; period, model_function = q_inverted_model);

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
instance_0 = zeroth_model(MeshParams(), model_0, channels_0);

# ╔═╡ 4dd7d6f5-0711-4156-b2e1-20a017d5276c
instance_q_uniform = q_uniform_model(MeshParams(), model_q_uniform, channels_0);

# ╔═╡ 351a6436-557d-4d76-871e-1f9e466de676
instance_q_inv = q_inverted_model(MeshParams(), model_q_inverted, channels_0);

# ╔═╡ f278a1e6-f459-4582-9400-9bfa2d852d14
begin
	local m_giant = 0.8
	local m_dwarf = 1.0
	local mass_quotient = m_dwarf / m_giant
	local mass_quotient_inv = m_giant / m_dwarf

	local cos_i = 0.5
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

	optim_result = maximize(logpdf_objective(model), x0)
	optim_array = Optim.maximizer(optim_result)

	var_info = DynamicPPL.unflatten(var_info, optim_array)
	tuples = Turing.Inference.getparams(model, var_info)
	return Chains(reshape(last.(tuples), (1, :, 1)), Symbol.(first.(tuples)))
end

# ╔═╡ f7c074e7-d5b1-4b2c-9ab6-57e85b237a04
MAP_0 = MAP_point(instance_0, init_params_0)

# ╔═╡ 67391953-6f49-466e-bc72-349d9df43b0b
MAP_q_uniform = MAP_point(instance_q_uniform, init_params_0)

# ╔═╡ 09c08a94-a567-4ba8-af64-424e9dbbe75b
MAP_q_inv = MAP_point(instance_q_inv, init_params_0)

# ╔═╡ d1e094f8-a0fd-45ce-9255-4b9de702433b
begin
	p = plot_points_days(channels_0, period)
	plot_lines_days(model_q_inverted, interpolated_mesh, channels_0, [MAP_q_inv], p; label = "q ~ Uniform")
	plot_lines_days(model_q_uniform, interpolated_mesh, channels_0, [MAP_q_uniform], p; label = "q^-1 ~ Uniform")
	plot_lines_days(model_0, interpolated_mesh, channels_0, [MAP_0], p; label = "m_giant, m_dwarf ~ Uniform")
	plot!(legend = true)
end

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

# ╔═╡ 130e2666-8c6b-4c9a-bb6b-00a77cd6ac01
WaldWolfowitzTest(differences(model_q_inverted, channels_0, MAP_q_inv)[1])

# ╔═╡ 6ebf9f3b-f877-45b6-af9d-715ef5b88b4e
WaldWolfowitzTest(differences(model_q_uniform, channels_0, MAP_q_uniform)[1])

# ╔═╡ Cell order:
# ╠═11decfb6-8f0a-11ee-1bf8-d3faf8756b8b
# ╠═6b6cf7e6-340f-4e34-aeef-c7c9f1695380
# ╠═0a209f5e-bd22-4072-9247-c0ddcd42fdb8
# ╠═18694e43-32c0-424f-ade0-63c6fc8af923
# ╠═4b30de9a-e872-422c-9b12-2cc67942aba3
# ╠═49216477-bda7-487c-9c85-56e1ed842d84
# ╠═c73161c8-0eb3-40af-b527-40a05e10ae24
# ╠═67b9b67f-7b71-4e16-904f-c89e244bea7e
# ╠═a6dccb45-fa0c-4aff-ae2e-23efc08fd8da
# ╠═e36d9949-98cc-4635-835c-400a584fcde4
# ╠═a10f8c4b-7d02-43d4-b5cf-e869e794c3d6
# ╠═471730f6-cd54-45cb-871d-dfdc7d7e4d65
# ╠═16711745-a609-4add-8118-3918af8293aa
# ╠═7773cf80-05cc-4c4d-9185-a81b1a237fa1
# ╠═15fa2682-e8e7-44b2-993a-3d1c7f66cf4e
# ╠═fca42ad4-ef24-4480-8573-6f5b3dab8584
# ╠═d2bf9a47-5c7a-4db8-b028-8d3f7383ff95
# ╠═4dd7d6f5-0711-4156-b2e1-20a017d5276c
# ╠═351a6436-557d-4d76-871e-1f9e466de676
# ╠═5449571a-5e90-4686-84de-6bfe09fba48f
# ╠═f278a1e6-f459-4582-9400-9bfa2d852d14
# ╠═535fb7b7-781c-4a48-b5f4-64aa81ff8bed
# ╠═b1701d19-4cb7-4c35-943e-6a9cd415c4a1
# ╠═a3b38672-7e73-4701-abc6-30d057f09754
# ╠═f7c074e7-d5b1-4b2c-9ab6-57e85b237a04
# ╠═67391953-6f49-466e-bc72-349d9df43b0b
# ╠═09c08a94-a567-4ba8-af64-424e9dbbe75b
# ╠═d1e094f8-a0fd-45ce-9255-4b9de702433b
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
# ╠═130e2666-8c6b-4c9a-bb6b-00a77cd6ac01
# ╠═6ebf9f3b-f877-45b6-af9d-715ef5b88b4e
