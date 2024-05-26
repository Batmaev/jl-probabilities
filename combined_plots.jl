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
model_q_uniform = ModelParams(; period, model_function = q_uniform_model);

# ╔═╡ 1baeb4e6-184b-4ae6-8fb1-f381651c2f39
model_q_inverted = ModelParams(; period, model_function = q_inverted_model);

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

# ╔═╡ 5fcb704e-eb42-42ea-825b-2f3a572f5c75
samples_0 = load_cache("0edb45400843111e1b569d0380ffefb79a85df60");

# ╔═╡ d1023bb5-0b12-470c-85ca-43a728751671
samples_q_uniform = load_cache("bb223212417baaf6ff519930e524a6fb35446450");

# ╔═╡ c644dbd3-3781-40ba-bd25-93a44a7d3771
samples_q_inv = load_cache("941b99e8d4a7c02a9ac23a04bbd56620e468f59b");

# ╔═╡ 8ca30dfd-18df-4048-bd60-0ceffd802648
begin
	gr()
	mass_ratio_plot = density(title = L"q = m_{\mathrm{giant}} / m_{\mathrm{dwarf}}")
	density!(samples_q_inv[:mass_quotient_inv], label = L"q\ \ \ \ \sim \mathrm{Uniform}")
	density!(samples_q_uniform[:mass_quotient_inv], label = L"q^{-1} \sim \mathrm{Uniform}")
	density!(samples_0[:mass_quotient_inv], label = L"m_{\mathrm{giant}}, m_{\mathrm{dwarf}} \sim \mathrm{Uniform}")
end

# ╔═╡ 423026e4-2cc6-461d-b225-9997a12019bf
savefig(mass_ratio_plot, "tex/pic_drafts/mass_ratio.pdf")

# ╔═╡ d48474e5-d8fe-4a03-9fde-87e5f061a362
begin
	gr()
	inclination_plot = density(title = "Наклонение (°)")
	density!(rad2deg.(samples_q_inv[:observer_angle]), label = L"q\ \ \ \ \sim \mathrm{Uniform}")
	density!(rad2deg.(samples_q_uniform[:observer_angle]), label = L"q^{-1} \sim \mathrm{Uniform}")
	density!(rad2deg.(samples_0[:observer_angle]), label = L"m_{\mathrm{giant}}, m_{\mathrm{dwarf}} \sim \mathrm{Uniform}")
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
			markersize = 2,
		)
	end

	scatter!(; kwargs...)
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
			magnitudes = star_magnitude(phases, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]

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
			markersize = 2;
			kwargs...
		)
		scatter!(
			subplot,
			phases .- 1,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2,
		)
		scatter!(
			subplot,
			phases .+ 1,
			channel.measurements_y,
			yerr = channel.σ_measured,
			markercolor = 1,
			markersize = 2,
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

# ╔═╡ 592a603c-67c5-4b7e-8fda-98f8ba484e06
plotlyjs()

# ╔═╡ 4d603280-5540-4d5f-8e36-0a59067d22f1
begin
	local p = plot_ribbon_phases(model_0, mesh, channels, samples_0; color = 1, label = "m_giant, m_dwarf ~ Uniform")
	plot_ribbon_phases(model_q_uniform, mesh, channels, samples_q_uniform, p; color = 2, label = "q^-1 ~ Uniform")
	plot_ribbon_phases(model_q_inverted, mesh, channels, samples_q_inv, p; color = 3, label = "q ~ Uniform")
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

# ╔═╡ 0d3ba047-c50b-40af-9aaf-bf2a0b4a8091
MAP_q_uniform = get_optim_params(samples_q_uniform)

# ╔═╡ fdca93ca-503c-4bd2-9a05-5777670221a1
MAP_q_inv = get_optim_params(samples_q_inv)

# ╔═╡ 530ea0d2-0965-4826-88e2-2bd5da5a0956
begin
	local p = plot_points_days(channels, period)
	plot_lines_days(model_0, mesh, channels, [MAP_0], p; label = "m_giant, m_dwarf ~ Uniform")
	plot_lines_days(model_q_uniform, mesh, channels, [MAP_q_uniform], p; label = "q^-1 ~ Uniform")
	plot_lines_days(model_q_inverted, mesh, channels, [MAP_q_inv], p; label = "q^-1 ~ Uniform")
	#plot!(legend = :outerright)
end

# ╔═╡ Cell order:
# ╠═4a348134-1a98-11ef-290d-933d7a2a61ac
# ╠═b1d9141a-c674-4087-aaa6-d372fe0c8270
# ╠═b8772986-d938-4bc4-b71a-4d0a393e8032
# ╠═92329f85-a715-4453-b4da-910b7c6786bc
# ╠═ec0ba187-904b-4253-b7d9-82815cd4af58
# ╠═1b51f03d-fdef-4a1d-aa6a-2d0a352db9f9
# ╠═1baeb4e6-184b-4ae6-8fb1-f381651c2f39
# ╠═edf1a809-72e9-4b81-9e73-2abf3658210d
# ╠═5fcb704e-eb42-42ea-825b-2f3a572f5c75
# ╠═d1023bb5-0b12-470c-85ca-43a728751671
# ╠═c644dbd3-3781-40ba-bd25-93a44a7d3771
# ╠═8ca30dfd-18df-4048-bd60-0ceffd802648
# ╠═423026e4-2cc6-461d-b225-9997a12019bf
# ╠═d48474e5-d8fe-4a03-9fde-87e5f061a362
# ╠═37e51e62-9dbc-4c8d-baf4-30530f446b88
# ╠═eb496939-cadc-4b6b-abfa-7322be698553
# ╠═27597a9f-56a9-45d4-af43-233fbc781638
# ╠═043df726-a3fa-47b9-85a9-1baf1809aaed
# ╠═d196a68f-385b-4da4-a624-cd520fe19b97
# ╠═acb8a0f7-c451-4b08-b8c0-a8ef7589890a
# ╠═592a603c-67c5-4b7e-8fda-98f8ba484e06
# ╠═4d603280-5540-4d5f-8e36-0a59067d22f1
# ╠═ea82eae9-e224-4598-9d58-001772de663d
# ╠═43d87082-c09e-4ff6-8c5c-33a742f7d5ba
# ╠═b9fa8af8-4803-4438-84f0-6e01bbf6873b
# ╠═0d3ba047-c50b-40af-9aaf-bf2a0b4a8091
# ╠═fdca93ca-503c-4bd2-9a05-5777670221a1
# ╠═530ea0d2-0965-4826-88e2-2bd5da5a0956
