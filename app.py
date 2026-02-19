import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- CORE SIMULATION ENGINE ---
class AbsorberSimulation:
    def __init__(self, gas_flow, gas_temp, target_conc, cw_temp):
        self.gas_flow = gas_flow  # kg/h
        self.gas_temp = gas_temp  # C
        self.hcl_frac = 0.6728    # Mass fraction in gas (20% Cl2, 80% HCl by mole)
        self.cw_temp = cw_temp
        self.target_conc = target_conc
        
        # Graphite Tube Geometry
        self.Di = 0.020  # ID meters
        
    def get_equilibrium_pressure(self, temp, conc):
        if conc < 0.01: return 0
        log_p = 10.5 - (2500 / (temp + 273)) + 2 * np.log(conc)
        return np.exp(log_p) * 133.322 # Convert to Pa

    def solve(self):
        total_hcl_in = self.gas_flow * self.hcl_frac
        product_rate = total_hcl_in / self.target_conc
        water_rate = product_rate * (1 - self.target_conc)
        
        # Hydraulics & Tube Count
        rho_gas = 2.42 
        u_gas_design = 16.0 # m/s
        vol_gas = (self.gas_flow / 3600) / rho_gas
        self.N_tubes = int(np.ceil((vol_gas / u_gas_design) / (np.pi * (self.Di/2)**2)))
        
        # Finite Difference Integration
        steps, L_step = 60, 0.1
        z_axis, temp_liquid, conc_liquid, flux_profile = [], [], [], []
        
        T_curr, C_curr = self.gas_temp, 0.05 
        current_hcl_mass = (water_rate * 0.05) / 0.95
        current_water_mass = water_rate
        accumulated_length, total_heat_removed = 0, 0

        for _ in range(steps):
            gamma = (current_water_mass + current_hcl_mass) / 3600 / (np.pi * self.Di * self.N_tubes)
            Re_L = 4 * gamma / 0.0008
            h_film = max(500, 1200 * (Re_L**(-1/3)))
            
            P_gas_partial = (current_hcl_mass / (current_hcl_mass + (self.gas_flow/3600))) * 101325 
            P_star = self.get_equilibrium_pressure(T_curr, C_curr)
            driving_force = max(0, P_gas_partial - P_star)
            
            mass_absorbed = (2.5e-4 * driving_force) * (np.pi * self.Di * L_step * self.N_tubes) * 36.5 / 1000 
            Q_gen = mass_absorbed * 1800 
            
            U_overall = 1.0 / (1/h_film + 1/2000 + 0.005/150) * 1000 
            Q_removed = U_overall * (np.pi * self.Di * L_step * self.N_tubes) * (T_curr - self.cw_temp) / 1000 
            
            dt = (Q_gen - Q_removed) / (((current_water_mass + current_hcl_mass)/3600) * 4.18)
            T_curr += dt
            current_hcl_mass += (mass_absorbed * 3600)
            C_curr = current_hcl_mass / (current_hcl_mass + current_water_mass)
            accumulated_length += L_step
            total_heat_removed += Q_removed
            
            z_axis.append(accumulated_length)
            temp_liquid.append(T_curr)
            conc_liquid.append(C_curr * 100)
            flux_profile.append(Q_removed)
            
            if C_curr >= self.target_conc: break
            
        return self.N_tubes, accumulated_length, total_heat_removed, max(temp_liquid), z_axis, temp_liquid, conc_liquid, flux_profile

# --- UI & DASHBOARD ---
st.set_page_config(page_title="HCl Absorber Simulation", layout="wide")

st.title("🧪 Advanced Falling Film Absorber Simulation")
st.markdown("Differential modeling of heat & mass transfer for HCl absorption in graphite tubes.")

# Sidebar Inputs
st.sidebar.header("Process Parameters")
gas_flow = st.sidebar.number_input("Gas Flow Rate (kg/h)", value=2500, step=100)
gas_temp = st.sidebar.number_input("Gas Inlet Temp (°C)", value=10, step=1)
target_conc = st.sidebar.number_input("Target HCl Conc (Mass Fraction)", value=0.32, step=0.01)
cw_temp = st.sidebar.number_input("Cooling Water Inlet Temp (°C)", value=32, step=1)

if st.sidebar.button("Run Simulation", type="primary"):
    with st.spinner("Solving differential equations..."):
        # Run calculation
        sim = AbsorberSimulation(gas_flow, gas_temp, target_conc, cw_temp)
        tubes, length, duty, max_t, z, t_liq, c_liq, flux = sim.solve()
        
        # Display KPIs
        st.subheader("Key Design Parameters")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Required Length", f"{length:.2f} m")
        col2.metric("Number of Tubes", f"{tubes}")
        col3.metric("Peak Temperature", f"{max_t:.2f} °C")
        col4.metric("Total Heat Duty", f"{duty:.2f} kW")
        
        # Plotly Charts
        st.subheader("Axial Profiles")
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=z, y=t_liq, name="Temperature (°C)", line=dict(color='firebrick', width=3)))
        fig.add_trace(go.Scatter(x=z, y=c_liq, name="Concentration (%)", yaxis='y2', line=dict(color='royalblue', width=3, dash='dash')))
        
        fig.update_layout(
            xaxis_title="Tube Length (m)",
            yaxis=dict(title="Temperature (°C)", titlefont=dict(color="firebrick"), tickfont=dict(color="firebrick")),
            yaxis2=dict(title="Concentration (%)", titlefont=dict(color="royalblue"), tickfont=dict(color="royalblue"), anchor="x", overlaying="y", side="right"),
            plot_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=40, r=40, t=40, b=40)
        )
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Heat Flux Profile")
        fig2 = go.Figure(data=[go.Bar(x=z, y=flux, marker_color='seagreen')])
        fig2.update_layout(xaxis_title="Tube Length (m)", yaxis_title="Heat Removed per segment (kW)", plot_bgcolor='rgba(0,0,0,0)')
        st.plotly_chart(fig2, use_container_width=True)
else:
    st.info("👈 Enter parameters in the sidebar and click 'Run Simulation' to see results.")
