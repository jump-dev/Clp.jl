# Copyright (c) 2013: Clp.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Clang.Generators
using Clp_jll

cd(@__DIR__)

clp_include_dir = joinpath(Clp_jll.artifact_dir, "include") |> normpath

coin_include_dir = joinpath(Clp_jll.CoinUtils_jll.artifact_dir, "include", "coin") |>normpath

options = load_options(joinpath(@__DIR__, "generate.toml"))

args = get_default_args()
push!(args, "-I$clp_include_dir", "-I$coin_include_dir")

headers = [
    joinpath(clp_include_dir, "coin", "Clp_C_Interface.h")
    joinpath(coin_include_dir, "Coin_C_defines.h")
]

ctx = create_context(headers, args, options)


build!(ctx)