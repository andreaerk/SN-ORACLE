print("Including Graphical Packages...")
using GeometryTypes
using GeometryBasics
using CairoMakie
using LaTeXStrings

print("Setting LaTeX Theme...")
set_theme!(theme_latexfonts())



function plot_settings(subplot; xlabel="X-label", ylabel="Y-label", 
                       xticks=nothing,  yticks=nothing,
                       xmticks=nothing, ymticks=nothing,
                       xscale = identity, yscale = identity,
                       title = "Title", titlesize = 10,
                       mirror_x = true, mirror_y = true,
                       xaxisposition=:bottom,
                       yaxisposition=:left,
                       xticklabelsize = 15, yticklabelsize = 15,
                       xticksize = 5, yticksize = 5,
                       xlabelsize=20, ylabelsize=20,
                       xticklabelsvisible=true,
                       yticklabelsvisible =true,
                       xlabelvisible=true, 
                       ylabelvisible=true,
                       xticklabelrotation=0, yticklabelrotation=0,
                       yminorticksize = 3,
                       xminorticksize = 3, 
                       yticklabelpad = 0,
                       ytrimspine = false,
                       ylabelpadding = 5,
                       xlabelpadding=5, 
                       xticklabelspace=Makie.automatic,
                       yticklabelspace=Makie.automatic,
       ) 

        

        if yscale== log10
                show_minor_y = true
                if yticks===nothing
                        yticks = 10. .^ collect(range(-99. ,+99. ,step=1))
                        yticks2 = ["$(@sprintf "%.0f" log10.(tick))" for tick in yticks]
                        yticks_t = [L"10^{%$(tick)}" for tick in yticks2]
                        yticks=(yticks, yticks_t)
        
                end 

                if ymticks ===nothing 
                        ymticks =  IntervalsBetween(9)
                end
        end

        if xscale== log10
                show_minor_x = true
                if xticks===nothing
                        xticks = 10. .^ collect(range(-99.,+99.,step=1))
                        xticks2 = ["$(@sprintf "%.0f" log10.(tick))" for tick in xticks]
                        xticks_t = [L"10^{%$(tick)}" for tick in xticks2]
                        xticks=(xticks, xticks_t)
        
                end 

                if xmticks===nothing
                        xmticks =  IntervalsBetween(9)
                end 
        end
        
        xticks===nothing ? xticks=Makie.automatic : xticks
        yticks===nothing ? yticks=Makie.automatic : yticks

        ytickformat=Makie.automatic
        xtickformat=Makie.automatic
        if xmticks!==nothing
                show_minor_x = true
        else
                show_minor_x = false 
                xmticks = [NaN]
        end
        if ymticks!==nothing
                show_minor_y = true
        else
                show_minor_y = false 
                ymticks = [NaN]
        end
        
        
        #@printf("yticks  %s \n", typeof(yticks))
        #@printf("ymticks %s \n", typeof(ymticks))
        #print(ymticks)
        #@printf("\n")
        #print(xticksize)
        #print(yticksize)

        ax = Axis(subplot, 
                xlabel = xlabel,  ylabel = ylabel, 
                title = title, titlesize = titlesize,
                xgridvisible = false, ygridvisible = false,
                xticks=xticks, yticks=yticks,
                xminorticks = xmticks, xminorticksvisible = show_minor_x,
                yminorticks = ymticks, yminorticksvisible = show_minor_y, 
                yminorticksize = yminorticksize, xminorticksize = xminorticksize, 
                yaxisposition = yaxisposition, xaxisposition = xaxisposition,
                xticksmirrored=mirror_x, yticksmirrored=mirror_y, 
                xtickalign=1, ytickalign=1, 
                xminortickalign=1, yminortickalign=1,
                xscale=xscale, yscale=yscale, 
                xtickformat=xtickformat, ytickformat=ytickformat,
                xticklabelsize=xticklabelsize, yticklabelsize=yticklabelsize,
                xlabelsize=xlabelsize, ylabelsize=ylabelsize,
                xticklabelsvisible=xticklabelsvisible,
                yticklabelsvisible=yticklabelsvisible,
                ylabelvisible=ylabelvisible, xlabelvisible=xlabelvisible,
                yticksize=yticksize, xticksize=xticksize,
                xticklabelrotation=xticklabelrotation, 
                yticklabelrotation=yticklabelrotation,
                yticklabelpad=yticklabelpad,xticklabelspace=xticklabelspace, yticklabelspace=yticklabelspace,
                ytrimspine = ytrimspine, ylabelpadding=ylabelpadding, xlabelpadding=xlabelpadding)

        return ax
end 


function plot_adjust(axis; xlabel="X-label", ylabel="Y-label", 
        xticks=nothing,  yticks=nothing,
        xmticks=nothing, ymticks=nothing,
        xscale = identity, yscale = identity,
        title = "Title", titlesize = 10,
        xaxisposition=:bottom,
        yaxisposition=:left,
        xticklabelsize = 15, yticklabelsize = 15,
        xticksize = 5, yticksize = 5,
        xlabelsize=20, ylabelsize=20,
        xticklabelsvisible=true,
        yticklabelsvisible =true,
        xlabelvisible=true,
        ylabelvisible=true,
        xticklabelrotation=0, yticklabelrotation=0,
        yminorticksize = 3,
        xminorticksize = 3, 
        rightspinevisible=true,
        leftspinevisible=true,
        topspinevisible=true,
        bottomspinevisible=true,
        xtickvisible = true,
        ytickvisible = true
        ) 

        ax.rightspinevisible = rightspinevisible
        ax.topspinevisible = topspinevisible
        ax.leftspinevisible = leftspinevisible
        ax.bottomspinevisible = bottomspinevisible
        ax.xticklabelsvisible = xticklabelsvisible
        ax.yticklabelsvisible = yticklabelsvisible
        ax.xticksvisible = xtickvisible
        ax.yticksvisible = ytickvisible
        ax.ylabelvisible = ylabelvisible
        ax.xlabelvisible = xlabelvisible
        return 0
end 

function cbar_settings(subplot; label="cbar-label", size=25,
        ticks=nothing, mticks=nothing, scale = identity,
        mirror = true,  axisposition=:right,
        ticksize = 15,labelsize=20, 
        colormap=cmap, colorrange=crange, highclip = :black, lowclip=:gray, ) 

        if scale== log10 && ticks===nothing
                ticks = 10. .^ collect(range(-99.,+99.,step=1))
                ticks2 = ["$(@sprintf "%.0f" log10.(tick))" for tick in ticks]
                ticks_t = [L"10^{%$(tick)}" for tick in ticks2]
                ticks=(ticks, ticks_t)
                mticks =  tick_logger(1e-99, 1e99)
                tickformat=Makie.automatic
                ticksize*=1.5
        else
                tickformat=Makie.automatic
        end

        ticks===nothing ? ticks=Makie.automatic : ticks

        if mticks!==nothing
                show_minor = true
        else
                show_minor = false 
                mticks = [NaN]
        end


        cbar = Colorbar(subplot, colormap=colormap, colorrange=colorrange, scale=scale, 
                        size=size, label=label, ticks=ticks,
                        highclip = highclip, lowclip = lowclip,
                        tickalign=1, labelsize=20, ticklabelsize=20, minorticks=mticks,
                        minorticksvisible=show_minor, minortickalign=1)

        return cbar
end 


"""
    make_hatch_pattern(; size, spacing, thickness, style, line_color, bg_color)

Returns a Matrix{RGBA{Float32}} with a hatch pattern drawn manually.

- `size`: image width and height (in pixels, can be Float64)
- `spacing`: distance between lines or dots
- `thickness`: width of lines or dot radius
- `style`: one of:
    - `:backslash` (`\`)
    - `:slash` (`/`)
    - `:cross` (`\` + `/`)
    - `:vertical`
    - `:horizontal`
    - `:plus` (`|` + `–`)
    - `:dotted`
- `line_color`: color of lines or dots
- `bg_color`: background color
"""
function make_hatch_pattern(; size=41, spacing=5, thickness=1.8,
                            style=:cross, line_color=RGBAf(0,0,0,1.0),
                            bg_color=RGBAf(1,1,1,0.0))

    w = Int(ceil(size))
    img = fill(bg_color, w, w)

    function draw_line(x0, y0, x1, y1)
        n = max(abs(x1 - x0), abs(y1 - y0)) * 2
        for i in range(0, 1, length=Int(n))
            x = round(Int, x0 + i*(x1 - x0))
            y = round(Int, y0 + i*(y1 - y0))
            for dx in -floor(Int, thickness/2):ceil(Int, thickness/2)
                for dy in -floor(Int, thickness/2):ceil(Int, thickness/2)
                    xi, yi = x + dx, y + dy
                    if 1 <= xi <= w && 1 <= yi <= w
                        img[yi, xi] = line_color
                    end
                end
            end
        end
    end

    if style in (:backslash, :cross)
        for offset in -w:spacing:2w
            draw_line(0, offset, offset, 0)
        end
    end

    if style in (:slash, :cross)
        for offset in -w:spacing:2w
            draw_line(0, w - offset, offset, w)
        end
    end

    if style in (:horizontal, :plus)
        for y in 0:spacing:w
            draw_line(0, y, w, y)
        end
    end

    if style in (:vertical, :plus)
        for x in 0:spacing:w
            draw_line(x, 0, x, w)
        end
    end

    if style == :dotted
        for x in 0:spacing:w
            for y in 0:spacing:w
                for dx in -floor(Int, thickness/2):ceil(Int, thickness/2)
                    for dy in -floor(Int, thickness/2):ceil(Int, thickness/2)
                        xi, yi = x + dx, y + dy
                        if 1 <= xi <= w && 1 <= yi <= w
                            img[yi, xi] = line_color
                        end
                    end
                end
            end
        end
    end

    return img
end


function stacked_barplot!(ax, x::AbstractVector, ydata;
                          width=0.8, direction=:y, colors=[], kwargs...)

    n = length(x)

    # Convert ydata to list-of-arrays format
    y_list = if isa(ydata, AbstractMatrix)
        [view(ydata, i, :) for i in 1:size(ydata, 1)]
    elseif isa(ydata, AbstractVector) && all(yi -> length(yi) == n, ydata)
        ydata
    else
        error("`y` must be a Matrix of size (s, n) or Vector of Vectors each of length $n")
    end

    s = length(y_list)

    # Handle color cycling
    if isempty(colors)
        colors = Makie.wong_colors()[1:s]
    end

    plots = Makie.Poly[]
    for j in 1:n  # bins
        base = 0.0
        for i in 1:s  # stack levels
            yij = y_list[i][j]
            if yij == 0
                continue
            end

            xc = x[j]
            if direction == :y
                rect = GeometryBasics.Rect(xc - width/2, base, width, yij)
                verts = [
                    Makie.Point2f(rect.origin),
                    Makie.Point2f(rect.origin[1] + GeometryBasics.width(rect), rect.origin[2]),
                    Makie.Point2f(rect.origin[1] + GeometryBasics.width(rect), rect.origin[2] + GeometryBasics.height(rect)),
                    Makie.Point2f(rect.origin[1], rect.origin[2] + GeometryBasics.height(rect)),
                ]
            elseif direction == :x
                rect = GeometryBasics.Rect(base, xc - width/2, yij, width)
                verts = [
                    Makie.Point2f(rect.origin),
                    Makie.Point2f(rect.origin[1] + GeometryBasics.width(rect), rect.origin[2]),
                    Makie.Point2f(rect.origin[1] + GeometryBasics.width(rect), rect.origin[2] + GeometryBasics.height(rect)),
                    Makie.Point2f(rect.origin[1], rect.origin[2] + GeometryBasics.height(rect)),
                ]
            else
                error("Unsupported direction $direction; use :x or :y")
            end

            color_val = colors[i % length(colors) == 0 ? length(colors) : i % length(colors)]
            p = Makie.poly!(ax, verts, color=color_val, kwargs...)
            push!(plots, p)
            base += yij
        end
    end

    return plots
end



"""
    hatched_pie!(ax, values; colors, radius=0.4, inner_radius=0.0, offset=π/2)

Draws a pie chart using poly! calls for each slice, allowing heterogeneous colors (e.g., hatch patterns).

- `ax`: Makie Axis to draw on.
- `values`: vector of non-negative values (one per slice).
- `colors`: vector of colorants, hatch patterns, or textures (same length as `values`).
- `radius`: outer radius of pie.
- `inner_radius`: if > 0, makes a donut chart.
- `offset`: starting angle in radians (default π/2 for 12 o'clock).
"""
function hatched_pie!(ax, values::AbstractVector;
                      colors,
                      radius=0.4,
                      inner_radius=0.0,
                      offset=π/2,
                      strokewidth=0.0,
                      strokecolor=:black)

    if any(values .< 0)
        error("All values must be non-negative")
    end
    if length(values) != length(colors)
        error("Length of `values` and `colors` must match")
    end

    total = sum(values)
    θ0 = offset
    center = Point2f(0, 0)

    for idx in eachindex(values)
        v = values[idx]
        c = colors[idx]  # Extract individually to avoid conversions

        if v == 0
            continue
        end

        θ1 = θ0 + 2π * v / total
        n = 50
        θs = range(θ0, θ1, length=n)

        outer = [center + radius * Point2f(cos(θ), sin(θ)) for θ in θs]

        if inner_radius == 0
        # Triangle fan: center point followed by arc
                verts = [center]
                append!(verts, outer)
        else
                inner = [center + inner_radius * Point2f(cos(θ), sin(θ)) for θ in reverse(θs)]
                verts = vcat(outer, inner)
        end

        # Use poly! per slice with individual color — avoids list conversion
        poly!(ax, verts; color=c, strokewidth=strokewidth, strokecolor=strokecolor)

        θ0 = θ1
    end

    return nothing
end



function move_axislabel!(axis; which = :x, align = 1.0)
    0 <= align <= 1 || @warn "Halign outside of (0, 1)!  Be warned."

    label_pos=nothing
    if which == :x
        label_pos = @lift $(axis.xaxis.elements[:labeltext][1])[1][1] # get the label position in the x dimension
    elseif which == :y
        label_pos = @lift $(axis.yaxis.elements[:labeltext][1])[1][1] # get the label position in the x dimension
    else
       throw(ErrorException)  
    end

    lift(axis.scene.viewport, label_pos) do px_area, current_pos
        xwidth = px_area.widths[1] # the width of the scene
        desired_offset = xwidth * align
        which == :x && translate!(axis.xaxis.elements[:labeltext], desired_offset - current_pos, 0, 0)
        which == :y && translate!(axis.yaxis.elements[:labeltext], 0, desired_offset - current_pos, 0)
    end
end


print("done!\n")
