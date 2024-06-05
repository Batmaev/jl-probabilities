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

	using Plots
	using StatsPlots

	using KernelDensity
	import PyPlot
	using Roots
	using Optim
end

# ╔═╡ cf1f5fab-fc8c-4e51-b364-c2cd56e002c2
using StatsBase

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
]

# ╔═╡ 6a91ffe3-d906-4269-9684-ac8ffa4f52c4
begin
	local m_giant = 0.8
	local m_dwarf = 1.0
	local mass_quotient = m_dwarf / m_giant

	local cos_i = 0.5
	local observer_angle = acos(cos_i)

	local initial_phase = -1.46867
	local offset = [46.6894, 48.8942]
	local log_σ_common = log.([0.01, 0.01])

	init_params = (; m_giant, m_dwarf, mass_quotient, cos_i, observer_angle, initial_phase, offset, log_σ_common)
end

# ╔═╡ 8db3b68e-1546-4c28-b7ff-b0ea57684fa7
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
	n_chains=8,
	n_samples = 120*32
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

# ╔═╡ ee35822f-7417-4d48-b799-1751d1f76f8f
(samples.info.stop_time - samples.info.start_time) / length(samples)

# ╔═╡ de408e58-efb8-4c7e-9415-d6e38e747d3f
begin
	sampled_values = samples[collect(values(samples.info.varname_to_symbol))]
	CSV.write("samples/0.csv", sampled_values)
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
		vals = Array{Float64}(undef, size(samples)[1], size(samples)[3], length(phases))

		for s1 ∈ 1 : size(samples)[1], s2 ∈ 1 : size(samples)[3]
			sample = group_symbols(samples[s1, :, s2])

			vals[s1, s2, :] = star_magnitude(phases .* 2π, interpolated_mesh, model_params, channel, sample) .+ sample[:offset][c]
		end

		means = mean(vals, dims = (1, 2))[:]
		stds = std(vals, dims = (1, 2))[:]

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
	rplot = plot_ribbon_phases(model_params, mesh, channels, sample(samples, 5_000), p)
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

# ╔═╡ 7853f148-3cfc-4def-946b-eadac15159ad
ess(samples)

# ╔═╡ c0983d9b-5bb5-438a-942b-14a8582efbad
samples[2, :, 3]

# ╔═╡ bc02da53-f9aa-439d-884e-87cd39da7f1f
PyPlot.svg(true)

# ╔═╡ 91bca842-224b-42a2-9ad6-6b337366e474
function biplot(x, y, levels, color = "black", ax = nothing)
	if ax == nothing
		PyPlot.figure()
		ax = PyPlot.gca()
	end

	k = kde((
		reshape(x, :),
		reshape(y, :),
	))

	mx = maximize(x -> pdf(k, x[1], x[2]), [mean(x), mean(y)])
	max_point = Optim.maximizer(mx)

	thresholds = get_threshold.(Ref(k), levels)

	ax.contour(
		k.x, k.y, k.density',
		levels = thresholds,
		colors = color,
	)
	ax.scatter([max_point[1]], [max_point[2]], color = color)
	print(max_point)
	PyPlot.gcf()
end

# ╔═╡ 5158fa88-c07a-4406-9f30-7ac5814ed8a2
biplot(samples[:mass_quotient_inv][:], rad2deg.(samples[:observer_angle][:]), [0.95, 0.68, 0])

# ╔═╡ 94df2966-5278-4364-badd-f9c56315b259
mean(samples[:mass_quotient_inv]), std(samples[:mass_quotient_inv])

# ╔═╡ ae5ce8b1-ebcf-4f86-b6c7-bc5b8bd8e864
function estimate_MAP_ci(arr, confidence_level)
	dist = kde(arr)
	mx = maximize(x -> pdf(dist, x[1]), [mean(arr)])
	mx = Optim.maximizer(mx)[1]

	total = sum(dist.density)
	function percentage(threshold)
		return sum(filter(pdf -> pdf > threshold, dist.density)) / total
	end

	threshold = find_zero(
		threshold -> percentage(threshold) - confidence_level,
		extrema(dist.density),
	)

	left_idx = findfirst(dist.density .> threshold)
	right_idx = findlast(dist.density .> threshold)

	return (mx, (dist.x[left_idx], dist.x[right_idx]))
end

# ╔═╡ d29c9e6e-2a5c-435e-8895-01b1cb2ff31c
est_q = estimate_MAP_ci(samples[:mass_quotient_inv].data[:], 0.68)

# ╔═╡ 621253f6-68d2-4d9e-8fb5-f4513a3a455b
est_q[2][1] - est_q[1], est_q[2][2] - est_q[1]

# ╔═╡ 91a362e6-1097-4ee1-8a45-b26d5261b14d
est_i = estimate_MAP_ci(rad2deg.(samples[:observer_angle].data[:]), 0.68)

# ╔═╡ 9e78fc3c-7ba2-46af-88dd-a1d5527aad8e
est_i[2][1] - est_i[1], est_i[2][2] - est_i[1]

# ╔═╡ 015e2ad0-65e9-43d3-b639-e1d77a348a28
biplot(samples[:m_dwarf][:], samples[:m_giant][:], [0.95, 0.68, 0])

# ╔═╡ 01d21747-7701-4658-b7f8-26d285af7bf3
biplot(samples[:m_dwarf][:], samples[:m_giant][:], [0.95, 0.68, 0])

# ╔═╡ ec91cfe1-9422-4017-ae9b-553a7bc86fdf
0.95 / 1.31

# ╔═╡ 7bf7d31b-d098-4bdd-b3c7-768c1c6c48c8
mode(samples[:m_dwarf])

# ╔═╡ 325e3d23-ba4e-4d95-9dd7-c2b1e91b6db9
coeftable(model_params.model_function(mesh_params, model_params, channels))

# ╔═╡ 15b39d12-3ca0-4cf3-8387-c004f4ce787e
coeftable(samples)

# ╔═╡ 841359b6-f44d-447c-b4e0-5a3bc2cb9341
confint(StatisticalModel(samples))

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
# ╠═8db3b68e-1546-4c28-b7ff-b0ea57684fa7
# ╠═4d347ae9-9b88-4a87-9459-5fec942cbd39
# ╠═c330906f-36e8-4b53-b0fc-a4b632024169
# ╠═c9f2e1d9-4d96-4d96-bcbd-02c6778bc055
# ╠═20913587-a0d2-45e8-88b9-1b1ab532d45e
# ╠═60698009-5b82-44e8-b8ee-da6432d8a227
# ╠═4513a94d-e0af-4842-b9cb-03a3176acfe3
# ╠═0c11e705-9921-4343-8ca9-b54ed3499af2
# ╠═ee35822f-7417-4d48-b799-1751d1f76f8f
# ╠═de408e58-efb8-4c7e-9415-d6e38e747d3f
# ╠═994aeb01-a8fb-4c15-af3e-f367fb237ae8
# ╠═fa3c2c79-6de1-4c4e-b6ef-93983917779a
# ╠═3cbab8f9-9416-4a17-8cf8-d5873a67dc72
# ╠═f3735b8f-7d73-4bad-be3e-4012d5b14a10
# ╠═017a39e2-a149-4cc6-befa-438911b87ec6
# ╠═884089d1-a9ba-41d0-85a6-8ce7de83cb2f
# ╠═7853f148-3cfc-4def-946b-eadac15159ad
# ╠═c0983d9b-5bb5-438a-942b-14a8582efbad
# ╠═bc02da53-f9aa-439d-884e-87cd39da7f1f
# ╠═91bca842-224b-42a2-9ad6-6b337366e474
# ╠═5158fa88-c07a-4406-9f30-7ac5814ed8a2
# ╠═94df2966-5278-4364-badd-f9c56315b259
# ╠═ae5ce8b1-ebcf-4f86-b6c7-bc5b8bd8e864
# ╠═d29c9e6e-2a5c-435e-8895-01b1cb2ff31c
# ╠═621253f6-68d2-4d9e-8fb5-f4513a3a455b
# ╠═91a362e6-1097-4ee1-8a45-b26d5261b14d
# ╠═9e78fc3c-7ba2-46af-88dd-a1d5527aad8e
# ╠═015e2ad0-65e9-43d3-b639-e1d77a348a28
# ╠═01d21747-7701-4658-b7f8-26d285af7bf3
# ╠═ec91cfe1-9422-4017-ae9b-553a7bc86fdf
# ╠═7bf7d31b-d098-4bdd-b3c7-768c1c6c48c8
# ╠═cf1f5fab-fc8c-4e51-b364-c2cd56e002c2
# ╠═325e3d23-ba4e-4d95-9dd7-c2b1e91b6db9
# ╠═15b39d12-3ca0-4cf3-8387-c004f4ce787e
# ╠═841359b6-f44d-447c-b4e0-5a3bc2cb9341
