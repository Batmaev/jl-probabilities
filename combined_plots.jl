### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 4a348134-1a98-11ef-290d-933d7a2a61ac
begin
	import Pkg
	Pkg.activate(".")

	using Revise
	using TuringStars

	using DelimitedFiles
	using DataFrames

	using Turing

	using Plots
	using StatsPlots
	using LaTeXStrings

	using KernelDensity
	import PyPlot
	using Roots
	using Optim
end

# ╔═╡ cb98d1cd-13b0-48c5-ad46-5c67a87bea48
using Dates

# ╔═╡ e99507a6-5429-4655-bc0d-d0b117cb2584
using HypothesisTests

# ╔═╡ af5264a4-549e-4d24-b12a-55393c2feb49
using CSV

# ╔═╡ b1d9141a-c674-4087-aaa6-d372fe0c8270
begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

# ╔═╡ b8772986-d938-4bc4-b71a-4d0a393e8032
period = 227.5687

# ╔═╡ 92329f85-a715-4453-b4da-910b7c6786bc
mesh = InterpolatedRocheMesh(MeshParams());

# ╔═╡ ec0ba187-904b-4253-b7d9-82815cd4af58
model_0 = ModelParams(; period);

# ╔═╡ 1b51f03d-fdef-4a1d-aa6a-2d0a352db9f9
model_mq_uniform = ModelParams(; period, model_function = mq_uniform_model);

# ╔═╡ 1baeb4e6-184b-4ae6-8fb1-f381651c2f39
model_mq_inverted = ModelParams(; period, model_function = mq_inverted_model);

# ╔═╡ edf1a809-72e9-4b81-9e73-2abf3658210d
channels = [
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

# ╔═╡ b6abfda3-0006-4de3-b191-26a00999008f
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

# ╔═╡ 5fcb704e-eb42-42ea-825b-2f3a572f5c75
samples_0 = load_cache("cb62aa114352825ae90b30a7a86b92a48a5e5ebe");

# ╔═╡ d1023bb5-0b12-470c-85ca-43a728751671
samples_mq_uniform = load_cache("cc0e51c820dbe7a404ef069d63b80c637e285884");

# ╔═╡ c644dbd3-3781-40ba-bd25-93a44a7d3771
samples_mq_inv = load_cache("ddcf791d611ecfc100e8bd870bfbe2ac33bbfc40");

# ╔═╡ 4892c5ad-54ac-449a-a50c-b4146fb7eb91
samples_bb = load_cache("45c8116e1d6458b51cb4b16444e84467b76af300");

# ╔═╡ 8ca30dfd-18df-4048-bd60-0ceffd802648
begin
	gr()
	mass_ratio_plot = density(title = L"q = m_{\mathrm{giant}} / m_{\mathrm{dwarf}}")
	density!(samples_mq_inv[:mass_quotient_inv].data[:], label = L"q\ \ \ \ \sim \mathrm{Uniform}")
	density!(samples_mq_uniform[:mass_quotient_inv].data[:], label = L"q^{-1} \sim \mathrm{Uniform}")
	density!(samples_0[:mass_quotient_inv].data[:], label = L"m_{\mathrm{giant}}, m_{\mathrm{dwarf}} \sim \mathrm{Uniform}")
	density!(samples_bb[:mass_quotient_inv].data[:], label = L"blackbody; $q \sim \mathrm{Uniform}$")
end

# ╔═╡ 423026e4-2cc6-461d-b225-9997a12019bf
savefig(mass_ratio_plot, "tex/pic_drafts/mass_ratio.pdf")

# ╔═╡ d48474e5-d8fe-4a03-9fde-87e5f061a362
begin
	gr()
	inclination_plot = density(title = "Наклонение (°)")
	density!(rad2deg.(samples_mq_inv[:observer_angle].data[:]), label = L"q\ \ \ \ \sim \mathrm{Uniform}")
	density!(rad2deg.(samples_mq_uniform[:observer_angle].data[:]), label = L"q^{-1} \sim \mathrm{Uniform}")
	density!(rad2deg.(samples_0[:observer_angle].data[:]), label = L"m_{\mathrm{giant}}, m_{\mathrm{dwarf}} \sim \mathrm{Uniform}")
	density!(rad2deg.(samples_bb[:observer_angle].data[:]), label = L"blackbody; $q \mathrm{Uniform}$")
end

# ╔═╡ 37e51e62-9dbc-4c8d-baf4-30530f446b88
function group_symbols(sample)
	if isa(sample, Chains)
		sample = get_params(sample)
	end
	sample
end

# ╔═╡ eb496939-cadc-4b6b-abfa-7322be698553
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

# ╔═╡ 27597a9f-56a9-45d4-af43-233fbc781638
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
			markersize = 2;
			kwargs...
		)
	end
	p
end

# ╔═╡ 043df726-a3fa-47b9-85a9-1baf1809aaed
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
			magnitudes = star_magnitude_smoker(phases, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]

			plot!(
				subplot,
				days .% period,
				magnitudes;
				kwargs...
			)
		end
	end
	p
end

# ╔═╡ d196a68f-385b-4da4-a624-cd520fe19b97
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
			markersize = 2,
			label = false;
			kwargs...
		)
		scatter!(
			subplot,
			phases .- 1,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2,
			label = false;
			kwargs...
		)
		scatter!(
			subplot,
			phases .+ 1,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2,
			label = false;
			kwargs...
		)
	end
	scatter!(xlim = (-0.25, 1.25))
end

# ╔═╡ acb8a0f7-c451-4b08-b8c0-a8ef7589890a
function plot_ribbon_phases(model_params, interpolated_mesh,
							channels, samples, p = nothing; kwargs...)
	if p == nothing
		p = plot_template("Фаза")
	end

	phases = -0.25 : 0.01 : 1.25

	for (c, (channel, subplot)) ∈ enumerate(zip(channels, p.subplots))
		vals = Array{Float64}(undef, length(samples), length(phases))

		for s ∈ 1 : length(samples)
			sample = group_symbols(samples[s])

			vals[s, :] = star_magnitude(phases .* 2π, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]
		end

		means = mean(vals, dims = 1)[1, :]
		stds = std(vals, dims = 1)[1, :]

		plot!(
			subplot,
			phases,
			means,
			alpha = 0,
			ribbon = stds,
			color = 1;
			kwargs...
		)
	end

	return p
end

# ╔═╡ c9c6392a-54e7-458c-9123-b4a6a90bca0f
function plot_lines_phases(model_params, interpolated_mesh,
							channels, samples, p = nothing; kwargs...)
	if p == nothing
		p = plot_template("Фаза")
	end

	phases = -0.25 : 0.01 : 1.25

	for (c, (channel, subplot)) ∈ enumerate(zip(channels, p.subplots))

		for s ∈ 1 : length(samples)
			sample = group_symbols(samples[s])

			vals = star_magnitude_smoker(phases .* 2π, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]

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

# ╔═╡ 592a603c-67c5-4b7e-8fda-98f8ba484e06
plotlyjs()

# ╔═╡ 4d603280-5540-4d5f-8e36-0a59067d22f1
begin
	local p = plot_ribbon_phases(model_0, mesh, channels, sample(samples_0, 2000); color = 1, label = "m_giant, m_dwarf ~ Uniform")
	plot_ribbon_phases(model_mq_uniform, mesh, channels, sample(samples_mq_uniform, 2000), p; color = 2, label = "q^-1 ~ Uniform")
	plot_ribbon_phases(model_mq_inverted, mesh, channels, sample(samples_mq_inv, 2000), p; color = 3, label = "q ~ Uniform")
	plot!(legend = true)
end

# ╔═╡ ea82eae9-e224-4598-9d58-001772de663d
function estimate_MAP_1D(arr)
	dist = kde(arr)
	mx = maximize(x -> pdf(dist, x[1]), [mean(arr)])
	return Optim.maximizer(mx)[1]
end

# ╔═╡ 43d87082-c09e-4ff6-8c5c-33a742f7d5ba
function get_optim_params(samples)
	names = samples.name_map.parameters
	vals = [estimate_MAP_1D(samples[col].data[:]) for col ∈ names]
	return Chains(reshape(vals, (1, :, 1)), names)
end

# ╔═╡ b9fa8af8-4803-4438-84f0-6e01bbf6873b
MAP_0 = get_optim_params(samples_0)

# ╔═╡ e69c4ebe-5dcb-4a43-8147-b83218404db5
0.724759 / 1.37872

# ╔═╡ 0d3ba047-c50b-40af-9aaf-bf2a0b4a8091
MAP_mq_uniform = get_optim_params(samples_mq_uniform)

# ╔═╡ fdca93ca-503c-4bd2-9a05-5777670221a1
MAP_mq_inv = get_optim_params(samples_mq_inv)

# ╔═╡ 88bb2acf-f00c-430e-b06b-5c781b9c1dc8
rad2deg.((acos(MAP_mq_uniform[:cos_i][1]), MAP_mq_uniform[:observer_angle][1]))

# ╔═╡ ce74332a-56f6-488b-9cc1-51e2ed3c9268
1.43626 * 0.701986

# ╔═╡ be417621-4d7f-4ad5-b368-461c14fcbe7c
MAP_bb = get_optim_params(samples_bb)

# ╔═╡ 1bd81707-51bd-48d5-b3ea-7442e3c559b0
function fix_changed_variables(model_params, sample, channels = channels)
	model = model_params.model_function(MeshParams(), model_params, channels)
	var_info = DynamicPPL.VarInfo(model)

	param_array = sample.value.data[:]
	var_info = DynamicPPL.unflatten(var_info, param_array)

	tuples = Turing.Inference.getparams(model, var_info)
	return Chains(reshape(last.(tuples), (1, :, 1)), Symbol.(first.(tuples)))
end

# ╔═╡ d2c9ac7e-5bbd-4a35-8c89-a08ec5f6874c
fix_changed_variables(model_0, MAP_bb)

# ╔═╡ f2c1239a-509f-4cae-9ff9-d82f3df4dd1c


# ╔═╡ 59bb71c5-6acb-45cd-8933-d305b1eda32e
model_instance = model_0.model_function(MeshParams(), model_0, channels)

# ╔═╡ 855c2ff9-6682-4bb8-891b-c52055105558
var_info = DynamicPPL.VarInfo(model_instance)

# ╔═╡ 441e3647-16fd-4d1b-9b12-9aa43bf82b15
syms = DynamicPPL.syms(var_info)

# ╔═╡ 9c654065-723b-4433-bc3f-932d9e246d07
vcat([group_symbols(MAP_0)[sym] for sym ∈ syms]...)

# ╔═╡ 3fdb796b-65bb-44f4-b5b7-59804e4a9a2a
var_info

# ╔═╡ 714d586b-f065-4fab-a517-3aa0dec28155
DynamicPPL.syms(var_info)

# ╔═╡ ca435472-a6c5-4e21-9eb0-18e07541e93d
predict(model_instance, MAP_0)

# ╔═╡ 7982ab8a-4b81-43ca-a8d4-9407e22a3ae4
0.709109 / 1.34774

# ╔═╡ a416587a-6530-44ea-8e54-9e41a858121f
gr()

# ╔═╡ 530ea0d2-0965-4826-88e2-2bd5da5a0956
map_plot = begin
	local p = plot_points_phases(channels, period, mean(mean([MAP_mq_inv[:initial_phase], MAP_mq_uniform[:initial_phase], MAP_0[:initial_phase]])), markercolor = "black")

	plot_lines_phases(model_mq_inverted, mesh, channels, [MAP_mq_inv], p; label = L"q \sim U(0.1, 10)", linecolor = 1)

	plot_lines_phases(model_mq_uniform, mesh, channels, [MAP_mq_uniform], p; label = L"q^{-1} \sim U(0.1, 10)", linecolor = 2)

	plot_lines_phases(model_0, mesh, channels, [MAP_0], p; label = L"m_\mathrm{giant} \sim U(0.6 M_\odot, 10 M_\odot), m_\mathrm{dwarf} \sim U(0.3 M_\odot, 1.44 M_\odot)", linecolor = 3)
 
	plot!(plot_title = "MAP-оценки кривых блеска", plot_titlevspan = 0.065, bottom_margin = 0Plots.px)
	plot!(p.subplots[2], legend = :bottom, legend_title = "Априорное распределение q", legend_title_font_pointsize = 8)
end

# ╔═╡ 81fa3f6f-5baf-4d5b-a0e4-37e97b1a3b4f
savefig(map_plot, "tex/pic_drafts/map_different_q.svg")

# ╔═╡ b8a9ddce-6b91-49cd-a2c9-808aabcdd9d2
begin
	local p = plot_points_phases(channels, period, mean(mean([MAP_mq_inv[:initial_phase], MAP_mq_uniform[:initial_phase], MAP_0[:initial_phase]])), markercolor = "black")

	plot_lines_phases(model_mq_inverted, mesh, channels, [fix_changed_variables(model_mq_inverted, MAP_mq_inv)], p; label = L"q \sim U(0.1, 10)", linecolor = 1)

	plot_lines_phases(model_mq_uniform, mesh, channels, [fix_changed_variables(model_mq_uniform, MAP_mq_uniform)], p; label = L"q^{-1} \sim U(0.1, 10)", linecolor = 2)

	plot_lines_phases(model_0, mesh, channels, [fix_changed_variables(model_0, MAP_0)], p; label = L"m_\mathrm{giant} \sim U(0.6 M_\odot, 10 M_\odot), m_\mathrm{dwarf} \sim U(0.3 M_\odot, 1.44 M_\odot)", linecolor = 3)
 
	plot!(plot_title = "MAP-оценки кривых блеска", plot_titlevspan = 0.065, bottom_margin = 0Plots.px)
	plot!(p.subplots[2], legend = :bottom, legend_title = "Априорное распределение q", legend_title_font_pointsize = 8)
end

# ╔═╡ e3b69f0f-17de-4e3b-87ec-b4d2eaece7c0
md"### Runs test"

# ╔═╡ 5f4adfcf-fa63-4efd-8ca8-8f0fe353f5ef
function differences(model_params, channels, MAP_params)
	MAP_params = group_symbols(MAP_params)

	map(enumerate(channels)) do (i, channel)

		phases = channel.measurements_t ./ model_params.period .* 2π .+ MAP_params[:initial_phase]
		phases = reshape(phases, :)

		predicted_magnitudes = MAP_params[:offset][i]
		predicted_magnitudes = star_magnitude_smoker(phases, mesh, model_params, channel, MAP_params)
		predicted_magnitudes .+= MAP_params[:offset][i]

		@. channel.measurements_y - predicted_magnitudes
	end
end

# ╔═╡ 52178d30-0fa4-48a1-b3c4-be19e2bd6cfc
scatter(
	# julian2datetime.(2450199.5 .+ points.day),
	#(points.day .% estimated_period) / estimated_period,
	sign.(differences(model_0, channels, MAP_0)[1]),
	#label = ["K" "J"],
	legend = false,
	markersize = 2,
	size = (600, 200)
)

# ╔═╡ 5b56abe8-4543-447e-acef-6461336621b0
function plot_positive_negative(x, y; p = nothing, kwargs...)
	positive_mask = y .> 0
	negative_mask = .! positive_mask

	if p == nothing
		p = plot()
	end

	scatter!(x[negative_mask], y[negative_mask], markercolor = 1; kwargs...)
	scatter!(x[positive_mask], y[positive_mask], markercolor = 2; kwargs...)
end

# ╔═╡ 19150041-390a-480c-8c92-916eef117bd2
plot_positive_negative(y; p = nothing, kwargs...) = plot_positive_negative(1 : length(y), y; p, kwargs...)

# ╔═╡ 86234b48-c2df-4674-b0f5-eaf927fa5a8f
gr()

# ╔═╡ 35ec55c2-d6ad-4015-b65e-f57ce875db08
begin
	diffs_JK = differences(model_mq_inverted, channels, MAP_mq_inv)
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

# ╔═╡ 91e86026-83ee-4dd8-8e56-184feeece890
savefig(diffs_plot, "tex/pic_drafts/residuals.pdf")

# ╔═╡ 69f7f276-cc13-49d3-af78-b89b1e573689
scatter(
	# julian2datetime.(2450199.5 .+ points.day),
	# (points.day .% period) / period,
	sign.(differences(model_0, channels, MAP_0)[1]),
	#label = ["K" "J"],
	legend = false,
	markersize = 2,
	size = (600, 300),
	xlabel = "Номер точки",
	ylabel = L"m_\mathrm{exp} - m_\mathrm{model}",
)

# ╔═╡ 2fa84acc-09ce-439a-98a1-d36482d97c4a
.!(differences(model_0, channels, MAP_0)[1] .> 0)

# ╔═╡ 7b5b06cb-eca3-44e3-b4bd-db37f0e5e69b
WaldWolfowitzTest(differences(model_0, channels, MAP_0)[1])

# ╔═╡ a81a33e9-9d56-4b60-b77f-6ab5c34d980b
WaldWolfowitzTest(differences(model_mq_uniform, channels, MAP_mq_uniform)[1])

# ╔═╡ 1d3e2c60-4477-491c-b9ef-71f02cc63280
WaldWolfowitzTest(differences(model_mq_inverted, channels, MAP_mq_inv)[1])

# ╔═╡ 3508f691-c829-46f3-ac88-f1b82b6496ab
# begin
# 	sampled_values = samples_bb[collect(values(samples_bb.info.varname_to_symbol))]
# 	CSV.write("samples/black_body_0.csv", sampled_values)
# end;

# ╔═╡ Cell order:
# ╠═4a348134-1a98-11ef-290d-933d7a2a61ac
# ╠═b1d9141a-c674-4087-aaa6-d372fe0c8270
# ╠═b8772986-d938-4bc4-b71a-4d0a393e8032
# ╠═92329f85-a715-4453-b4da-910b7c6786bc
# ╠═ec0ba187-904b-4253-b7d9-82815cd4af58
# ╠═1b51f03d-fdef-4a1d-aa6a-2d0a352db9f9
# ╠═1baeb4e6-184b-4ae6-8fb1-f381651c2f39
# ╠═edf1a809-72e9-4b81-9e73-2abf3658210d
# ╠═b6abfda3-0006-4de3-b191-26a00999008f
# ╠═5fcb704e-eb42-42ea-825b-2f3a572f5c75
# ╠═d1023bb5-0b12-470c-85ca-43a728751671
# ╠═c644dbd3-3781-40ba-bd25-93a44a7d3771
# ╠═4892c5ad-54ac-449a-a50c-b4146fb7eb91
# ╠═8ca30dfd-18df-4048-bd60-0ceffd802648
# ╠═423026e4-2cc6-461d-b225-9997a12019bf
# ╠═d48474e5-d8fe-4a03-9fde-87e5f061a362
# ╠═37e51e62-9dbc-4c8d-baf4-30530f446b88
# ╠═eb496939-cadc-4b6b-abfa-7322be698553
# ╠═27597a9f-56a9-45d4-af43-233fbc781638
# ╠═043df726-a3fa-47b9-85a9-1baf1809aaed
# ╠═d196a68f-385b-4da4-a624-cd520fe19b97
# ╠═acb8a0f7-c451-4b08-b8c0-a8ef7589890a
# ╠═c9c6392a-54e7-458c-9123-b4a6a90bca0f
# ╠═592a603c-67c5-4b7e-8fda-98f8ba484e06
# ╠═4d603280-5540-4d5f-8e36-0a59067d22f1
# ╠═ea82eae9-e224-4598-9d58-001772de663d
# ╠═43d87082-c09e-4ff6-8c5c-33a742f7d5ba
# ╠═b9fa8af8-4803-4438-84f0-6e01bbf6873b
# ╠═e69c4ebe-5dcb-4a43-8147-b83218404db5
# ╠═0d3ba047-c50b-40af-9aaf-bf2a0b4a8091
# ╠═fdca93ca-503c-4bd2-9a05-5777670221a1
# ╠═88bb2acf-f00c-430e-b06b-5c781b9c1dc8
# ╠═ce74332a-56f6-488b-9cc1-51e2ed3c9268
# ╠═be417621-4d7f-4ad5-b368-461c14fcbe7c
# ╠═1bd81707-51bd-48d5-b3ea-7442e3c559b0
# ╠═d2c9ac7e-5bbd-4a35-8c89-a08ec5f6874c
# ╠═441e3647-16fd-4d1b-9b12-9aa43bf82b15
# ╠═9c654065-723b-4433-bc3f-932d9e246d07
# ╠═f2c1239a-509f-4cae-9ff9-d82f3df4dd1c
# ╠═59bb71c5-6acb-45cd-8933-d305b1eda32e
# ╠═855c2ff9-6682-4bb8-891b-c52055105558
# ╠═3fdb796b-65bb-44f4-b5b7-59804e4a9a2a
# ╠═714d586b-f065-4fab-a517-3aa0dec28155
# ╠═ca435472-a6c5-4e21-9eb0-18e07541e93d
# ╠═7982ab8a-4b81-43ca-a8d4-9407e22a3ae4
# ╠═a416587a-6530-44ea-8e54-9e41a858121f
# ╠═530ea0d2-0965-4826-88e2-2bd5da5a0956
# ╠═81fa3f6f-5baf-4d5b-a0e4-37e97b1a3b4f
# ╠═b8a9ddce-6b91-49cd-a2c9-808aabcdd9d2
# ╟─e3b69f0f-17de-4e3b-87ec-b4d2eaece7c0
# ╠═5f4adfcf-fa63-4efd-8ca8-8f0fe353f5ef
# ╠═52178d30-0fa4-48a1-b3c4-be19e2bd6cfc
# ╠═cb98d1cd-13b0-48c5-ad46-5c67a87bea48
# ╠═5b56abe8-4543-447e-acef-6461336621b0
# ╠═19150041-390a-480c-8c92-916eef117bd2
# ╠═86234b48-c2df-4674-b0f5-eaf927fa5a8f
# ╠═35ec55c2-d6ad-4015-b65e-f57ce875db08
# ╠═91e86026-83ee-4dd8-8e56-184feeece890
# ╠═69f7f276-cc13-49d3-af78-b89b1e573689
# ╠═2fa84acc-09ce-439a-98a1-d36482d97c4a
# ╠═e99507a6-5429-4655-bc0d-d0b117cb2584
# ╠═7b5b06cb-eca3-44e3-b4bd-db37f0e5e69b
# ╠═a81a33e9-9d56-4b60-b77f-6ab5c34d980b
# ╠═1d3e2c60-4477-491c-b9ef-71f02cc63280
# ╠═af5264a4-549e-4d24-b12a-55393c2feb49
# ╠═3508f691-c829-46f3-ac88-f1b82b6496ab
