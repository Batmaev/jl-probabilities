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
model_params = ModelParams(; period, model_function = mq_uniform_model)

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
	local mass_quotient = 1.38208

	local cos_i = 0.5
	local observer_angle = acos(cos_i)

	local initial_phase = -1.46867
	local offset = [46.6894, 48.8942]
	local log_σ_common = log.([0.01, 0.01])

	init_params = (; mass_quotient, cos_i, observer_angle, initial_phase, offset, log_σ_common)
end

# ╔═╡ 15055f84-78e5-416c-a9ad-a909f12b2f50
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
	n_chains = 8,
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

# ╔═╡ 97177b66-3f53-483e-a9e0-c31bf6852b5b
hash_chain_params(chain_params)

# ╔═╡ 0c11e705-9921-4343-8ca9-b54ed3499af2
samples = cached_sample(chain_params)

# ╔═╡ 2034d233-293a-4f5e-bcb3-5236a66e6699
(samples.info.stop_time - samples.info.start_time) / length(samples)

# ╔═╡ 802213be-600a-4f19-807d-bad83ac62257
begin
	sampled_values = samples[collect(values(samples.info.varname_to_symbol))]
	CSV.write("samples/mq_uniform.csv", sampled_values)
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
	plot_ribbon_phases(model_params, mesh, channels, sample(samples, 5_000), p)
end

# ╔═╡ 91bca842-224b-42a2-9ad6-6b337366e474
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

# ╔═╡ 7e14fcde-3895-4af5-8982-6623a62e2038
PyPlot.svg(true)

# ╔═╡ 2076e26c-12c1-4e34-a530-b13d71182f9e
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

# ╔═╡ c8f10134-2c3d-4271-baba-bc0abcb914cc
biplot(samples, [0.95, 0.68, 0])

# ╔═╡ 81907d0c-ed2d-4253-b252-9be3c841f3ee
mean(samples[:mass_quotient_inv]), std(samples[:mass_quotient_inv])

# ╔═╡ 6ac0c191-c563-47e9-a28f-b996dbb52364
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

# ╔═╡ c65d846c-8b09-4c2b-bbd7-4e91336f1029
est_q = estimate_MAP_ci(samples[:mass_quotient_inv].data[:], 0.68)

# ╔═╡ d4c95eb6-96ac-4e1f-b153-085747f4d8ec
est_q[2][1] - est_q[1], est_q[2][2] - est_q[1]

# ╔═╡ b5d6df83-80ae-46ed-819f-1deef9b8b25f
est_i = estimate_MAP_ci(rad2deg.(samples[:observer_angle].data[:]), 0.68)

# ╔═╡ 5412767a-c3bc-42c7-9755-1bfccb74a779
est_i[2][1] - est_i[1], est_i[2][2] - est_i[1]

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
# ╠═15055f84-78e5-416c-a9ad-a909f12b2f50
# ╠═4d347ae9-9b88-4a87-9459-5fec942cbd39
# ╠═c330906f-36e8-4b53-b0fc-a4b632024169
# ╠═c9f2e1d9-4d96-4d96-bcbd-02c6778bc055
# ╠═20913587-a0d2-45e8-88b9-1b1ab532d45e
# ╠═60698009-5b82-44e8-b8ee-da6432d8a227
# ╠═97177b66-3f53-483e-a9e0-c31bf6852b5b
# ╠═0c11e705-9921-4343-8ca9-b54ed3499af2
# ╠═2034d233-293a-4f5e-bcb3-5236a66e6699
# ╠═802213be-600a-4f19-807d-bad83ac62257
# ╠═994aeb01-a8fb-4c15-af3e-f367fb237ae8
# ╠═fa3c2c79-6de1-4c4e-b6ef-93983917779a
# ╠═3cbab8f9-9416-4a17-8cf8-d5873a67dc72
# ╠═f3735b8f-7d73-4bad-be3e-4012d5b14a10
# ╠═017a39e2-a149-4cc6-befa-438911b87ec6
# ╠═91bca842-224b-42a2-9ad6-6b337366e474
# ╠═7e14fcde-3895-4af5-8982-6623a62e2038
# ╠═2076e26c-12c1-4e34-a530-b13d71182f9e
# ╠═c8f10134-2c3d-4271-baba-bc0abcb914cc
# ╠═81907d0c-ed2d-4253-b252-9be3c841f3ee
# ╠═6ac0c191-c563-47e9-a28f-b996dbb52364
# ╠═c65d846c-8b09-4c2b-bbd7-4e91336f1029
# ╠═d4c95eb6-96ac-4e1f-b153-085747f4d8ec
# ╠═b5d6df83-80ae-46ed-819f-1deef9b8b25f
# ╠═5412767a-c3bc-42c7-9755-1bfccb74a779
