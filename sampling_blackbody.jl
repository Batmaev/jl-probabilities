### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 25ce51bc-19cb-11ef-22af-9d2f1925f669
begin
	import Pkg
	Pkg.activate(".")

	using Revise
	using TuringStars

	using DelimitedFiles
	using DataFrames
	using CSV

	using Turing

	using JSON3
	using SHA

	using Plots
	using StatsPlots

	using KernelDensity
	import PyPlot
	using Roots
	using Optim
end

# ╔═╡ 2c36d253-132b-471f-a395-28479e55562d
begin
	plotlyjs()
	theme(:juno)
end

# ╔═╡ f32fdc1c-272a-48c5-8f18-4a292d72f643
begin
	points = readdlm("stars/T_CrB_JK.dat")[2:end, :]
	points = map(x -> isa(x, Number) ? x : missing, points)
	points = DataFrame(points, [:day, :J, :J_err, :K, :K_err])
	points = dropmissing(points)
	points.day .-= points.day[1]
	points
end

# ╔═╡ 01b4866f-7ec3-4856-b7d8-adf8315da3ca
period = 227.5687

# ╔═╡ 93bda75e-fd00-4ac2-bd9f-a8045d4542c5
mesh_params = MeshParams()

# ╔═╡ 108f507b-e73e-474f-93ee-85cb9a006fc6
mesh = InterpolatedRocheMesh(mesh_params);

# ╔═╡ a6d96625-71e8-4b3a-b82e-9c21cdfdbab5
model_params = ModelParams(; period)

# ╔═╡ 0af70b03-c53f-415d-8f4f-f8bf6c383f3f
channels = [
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
]

# ╔═╡ 6a91ffe3-d906-4269-9684-ac8ffa4f52c4
begin
	local m_giant = 0.8
	local m_dwarf = 1.0
	local mass_quotient = m_dwarf / m_giant

	local cos_i = 0.5
	local observer_angle = acos(cos_i)

	local initial_phase = -1.46867
	local offset = [18.85, 21.17]
	local log_σ_common = log.([0.01, 0.01])

	init_params = (; m_giant, m_dwarf, mass_quotient, cos_i, observer_angle, initial_phase, offset, log_σ_common)
end

# ╔═╡ 81a1e646-8e76-44cb-8302-5ec00134599e
function group_symbols(sample)
	if isa(sample, Chains)
		sample = get_params(sample)
	end
	sample
end

# ╔═╡ c330906f-36e8-4b53-b0fc-a4b632024169
chain_params = ChainParams(;
	mesh_params,
	model_params,
	channels,
	init_params,
	n_samples = 4096
)

# ╔═╡ c9f2e1d9-4d96-4d96-bcbd-02c6778bc055
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

# ╔═╡ 20913587-a0d2-45e8-88b9-1b1ab532d45e
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

# ╔═╡ 60698009-5b82-44e8-b8ee-da6432d8a227
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
				magnitudes
			)
		end
	end

	plot!(; kwargs...)
end

# ╔═╡ 4d347ae9-9b88-4a87-9459-5fec942cbd39
begin
	local p = plot_points_days(channels, period)
	plot_lines_days(model_params, mesh, channels, [init_params], p)
end

# ╔═╡ 4513a94d-e0af-4842-b9cb-03a3176acfe3
hash_chain_params(chain_params)

# ╔═╡ 0c11e705-9921-4343-8ca9-b54ed3499af2
samples = cached_sample(chain_params)

# ╔═╡ 238b0807-35a4-40b9-b14e-1fabdd98cf3f
(samples.info.stop_time - samples.info.start_time) / length(samples)

# ╔═╡ e6821d90-c833-4b6c-9af1-433aba401571
begin
	sampled_values = samples[collect(values(samples.info.varname_to_symbol))]
	CSV.write("samples/black_body.csv", sampled_values)
end;

# ╔═╡ 994aeb01-a8fb-4c15-af3e-f367fb237ae8
begin
	local p = plot_points_days(channels, period)
	plot_lines_days(model_params, mesh, channels, sample(samples, 15), p)
end

# ╔═╡ fa3c2c79-6de1-4c4e-b6ef-93983917779a
plot(samples, margin = 10Plots.px, bottom_margin = 50Plots.px)

# ╔═╡ 3cbab8f9-9416-4a17-8cf8-d5873a67dc72
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
	scatter!(; kwargs...)
end

# ╔═╡ f3735b8f-7d73-4bad-be3e-4012d5b14a10
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
			color = 1
		)
	end

	plot!(; kwargs...)
end

# ╔═╡ 017a39e2-a149-4cc6-befa-438911b87ec6
begin
	local p = plot_points_phases(channels, period, mean(samples[:initial_phase]))
	plot_ribbon_phases(model_params, mesh, channels, samples, p)
end

# ╔═╡ 884089d1-a9ba-41d0-85a6-8ce7de83cb2f
function get_threshold(kde, confidence_level)
	total = sum(kde.density)
	function percentage(threshold)
		return sum(filter(pdf -> pdf > threshold, kde.density)) / total
	end
	return find_zero(
		threshold -> percentage(threshold) - confidence_level,
		extrema(kde.density),
	)
end

# ╔═╡ bc02da53-f9aa-439d-884e-87cd39da7f1f
PyPlot.svg(true)

# ╔═╡ 91bca842-224b-42a2-9ad6-6b337366e474
function biplot(samples, levels, color = "black", ax = nothing)
	if ax == nothing
		PyPlot.figure()
		ax = PyPlot.gca()
	end

	k = kde((
		reshape(1 ./ samples[:mass_quotient], :),
		reshape(rad2deg.(samples[:observer_angle]), :),
	))

	mx = maximize(x -> pdf(k, x[1], x[2]), [0.5, 50.])
	max_point = Optim.maximizer(mx)

	thresholds = get_threshold.(Ref(k), levels)

	ax.contour(
		k.x, k.y, k.density',
		levels = thresholds,
		colors = color,
	)
	ax.scatter([max_point[1]], [max_point[2]], color = color)
	PyPlot.gcf()
end

# ╔═╡ 5158fa88-c07a-4406-9f30-7ac5814ed8a2
biplot(samples, [0.95, 0.68, 0], "black")

# ╔═╡ Cell order:
# ╠═25ce51bc-19cb-11ef-22af-9d2f1925f669
# ╠═2c36d253-132b-471f-a395-28479e55562d
# ╠═f32fdc1c-272a-48c5-8f18-4a292d72f643
# ╠═01b4866f-7ec3-4856-b7d8-adf8315da3ca
# ╠═93bda75e-fd00-4ac2-bd9f-a8045d4542c5
# ╠═108f507b-e73e-474f-93ee-85cb9a006fc6
# ╠═a6d96625-71e8-4b3a-b82e-9c21cdfdbab5
# ╠═0af70b03-c53f-415d-8f4f-f8bf6c383f3f
# ╠═6a91ffe3-d906-4269-9684-ac8ffa4f52c4
# ╠═81a1e646-8e76-44cb-8302-5ec00134599e
# ╠═4d347ae9-9b88-4a87-9459-5fec942cbd39
# ╠═c330906f-36e8-4b53-b0fc-a4b632024169
# ╠═c9f2e1d9-4d96-4d96-bcbd-02c6778bc055
# ╠═20913587-a0d2-45e8-88b9-1b1ab532d45e
# ╠═60698009-5b82-44e8-b8ee-da6432d8a227
# ╠═4513a94d-e0af-4842-b9cb-03a3176acfe3
# ╠═0c11e705-9921-4343-8ca9-b54ed3499af2
# ╠═238b0807-35a4-40b9-b14e-1fabdd98cf3f
# ╠═e6821d90-c833-4b6c-9af1-433aba401571
# ╠═994aeb01-a8fb-4c15-af3e-f367fb237ae8
# ╠═fa3c2c79-6de1-4c4e-b6ef-93983917779a
# ╠═3cbab8f9-9416-4a17-8cf8-d5873a67dc72
# ╠═f3735b8f-7d73-4bad-be3e-4012d5b14a10
# ╠═017a39e2-a149-4cc6-befa-438911b87ec6
# ╠═884089d1-a9ba-41d0-85a6-8ce7de83cb2f
# ╠═bc02da53-f9aa-439d-884e-87cd39da7f1f
# ╠═91bca842-224b-42a2-9ad6-6b337366e474
# ╠═5158fa88-c07a-4406-9f30-7ac5814ed8a2
