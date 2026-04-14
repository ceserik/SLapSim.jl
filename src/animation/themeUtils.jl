function detect_dark_theme!()
    Sys.islinux() || return
    try
        xd = read(`kreadconfig6 --key LookAndFeelPackage`, String)
        is_dark = occursin("dark", lowercase(xd))
        println("Is dark theme: ", is_dark)
        if is_dark
            set_theme!(theme_dark())
        else
            set_theme!(Makie.current_default_theme())
        end
    catch
        println("Could not detect dark theme")
    end
end

function apply_slapsim_theme!()
    update_theme!(
        fontsize = 11,
        fonts = (; regular = "Latin Modern Roman", bold = "Latin Modern Roman Bold"),
        palette = (color = Makie.to_colormap(:tab10),),
        colormap = :turbo,
        Axis = (
            titlesize = 13,
            xlabelsize = 12,
            ylabelsize = 12,
            xticklabelsize = 10,
            yticklabelsize = 10,
        ),
        Legend = (
            labelsize = 11,
            titlesize = 12,
        ),
        Colorbar = (
            ticklabelsize = 10,
            labelsize = 11,
        ),
    )
end

function setup_plot_theme!()
    detect_dark_theme!()
    apply_slapsim_theme!()
end
