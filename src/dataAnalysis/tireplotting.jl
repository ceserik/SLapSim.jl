using SLapSim
using GLMakie

motor = createFischerMotor()
gearbox = createCTU25gearbox()

tire_lin = createR20lin(motor, gearbox)
tire_pac = SLapSim.createR20_pacejka(motor, gearbox)

Fz = 1000.0
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Slip angle [deg]", ylabel="Lateral force [N]", title="Tire comparison (Fz=$(Fz)N)")
SLapSim.plotTire(tire_lin; Fz=Fz, max_slip_deg=10.0, label="Linear", ax=ax)
SLapSim.plotTire(tire_pac; Fz=Fz, max_slip_deg=10.0, label="Pacejka", ax=ax)
axislegend(ax)
display(GLMakie.Screen(), fig)
