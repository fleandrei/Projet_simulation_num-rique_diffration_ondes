# code from https://genkuroki.github.io/documents/Jupyter/20170624%20Examples%20of%20animations%20in%20Julia%20by%20PyPlot%20and%20matplotlib.animation.html


# using PyCall
# @pyimport matplotlib.animation as anim
# using PyPlot
using Plots
using Makie

# function showanim(filename)
#     base64_video = base64encode(open(filename))
#     display("text/html", """<video controls src="data:video/x-m4v;base64,$base64_video">""")
# end

# A = randn(20,20,20)

# fig = figure(figsize=(2,2))

# function make_frame(i)
#     imshow(A[:,:,i+1],interpolation="none")
# end

# withfig(fig) do
#     myanim = anim.FuncAnimation(fig, make_frame, frames=size(A,3), interval=20)
#     myanim[:save]("test2.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
# end

# showanim("test2.mp4")

# 
# type(zeros(Float64, 5, 5))

## @gif or @animate macro should be followed by some sort of iteration
# A = randn(20,20,20)

# anim = @animate for i=1:10
#     p = rand(1, 10000)
#     histogram(p)
# end

# gif(anim, "mygif.gif", fps = 1)


# loadpath = "test"
# animation = Animation(loadpath,String[])
# u_sol = randn(20,20)
# for i=1:20
#         p = plot(u_sol[i,:],u_sol[i,:])
#         frame( animation, p )
# end
# name_of_gif = "test.mp4"
# run(`ffmpeg -framerate 15 -i $loadpath"%06d.png" $name_of_gif`) # run this in the REPL, it will hang indefinetly, if the gif already exists






 scene = lines(rand(10); linewidth=10)

 record(scene, "line_changing_colour.mp4", 1:255; framerate = 60) do i
        scene.plots[2][:color] = RGBf0(i/255, (255 - i)/255, 0) # animate scene
        # `scene.plots` gives the plots of the Scene.
        # `scene.plots[1]` is always the Axis if it exists,
        # and `scene.plots[2]` onward are the user-defined plots.
end