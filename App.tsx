import React, { useState, useMemo } from 'react';
import { 
  Rocket, 
  Settings, 
  Flame, 
  Ruler, 
  Cylinder, 
  Activity, 
  Play, 
  Box, 
  Droplets,
  Zap,
  Maximize,
  Gauge,
  CircleDashed,
  Wind,
  Thermometer,
  Weight,
  Timer,
  TrendingUp,
  Cone,
  ArrowUp,
  LayoutDashboard,
  LineChart as LineChartIcon,
  Lock,
  User,
  LogOut,
  Award,
  Mountain,
  Scale,
  Package
} from 'lucide-react';
import { 
  ComposedChart, 
  LineChart,
  Line, 
  XAxis, 
  YAxis, 
  CartesianGrid, 
  Tooltip, 
  ResponsiveContainer, 
  Area, 
  Legend, 
  AreaChart
} from 'recharts';

import { MotorInputs, ExposureOption, SimulationResults, SimulationStep, BurnRateStep } from './types';
import { INITIAL_INPUTS } from './constants';
import { InputGroup } from './components/InputGroup';
import { ResultCard } from './components/ResultCard';
import { MotorVisualizer } from './components/MotorVisualizer';

const getMotorClass = (impulse: number): string => {
  if (impulse < 1.26) return 'Micro';
  if (impulse < 2.5) return 'A';
  if (impulse < 5) return 'B';
  if (impulse < 10) return 'C';
  if (impulse < 20) return 'D';
  if (impulse < 40) return 'E';
  if (impulse < 80) return 'F';
  if (impulse < 160) return 'G';
  if (impulse < 320) return 'H';
  if (impulse < 640) return 'I';
  if (impulse < 1280) return 'J';
  if (impulse < 2560) return 'K';
  if (impulse < 5120) return 'L';
  if (impulse < 10240) return 'M';
  if (impulse < 20480) return 'N';
  if (impulse < 40960) return 'O';
  if (impulse < 81920) return 'P';
  if (impulse < 163840) return 'Q';
  if (impulse < 327680) return 'R';
  if (impulse < 655360) return 'S';
  if (impulse < 1310720) return 'T';
  if (impulse < 2621440) return 'U';
  if (impulse < 5242880) return 'V';
  return '>V';
};

const CLASS_ESTIMATES: Record<string, { altitude: string, weight: string, payload: string }> = {
  'A': { altitude: '0.1 - 0.7', weight: '0.0063 - 0.0125', payload: '0.000252 - 0.0005' },
  'B': { altitude: '0.4 - 2.9', weight: '0.0125 - 0.025', payload: '0.0005 - 0.001' },
  'C': { altitude: '1.7 - 11.6', weight: '0.025 - 0.05', payload: '0.001 - 0.002' },
  'D': { altitude: '7.0 - 46.5', weight: '0.05 - 0.1', payload: '0.002 - 0.004' },
  'E': { altitude: '27.9 - 186.1', weight: '0.1 - 0.2', payload: '0.004 - 0.008' },
  'F': { altitude: '111.7 - 744.4', weight: '0.2 - 0.4', payload: '0.008 - 0.016' },
  'G': { altitude: '446.7 - 2977.8', weight: '0.4 - 0.8', payload: '0.016 - 0.032' },
  'H': { altitude: '1786.6 - 11911', weight: '0.8 - 1.6', payload: '0.032 - 0.064' },
  'I': { altitude: '7146.8 - 47645', weight: '1.6 - 3.2', payload: '0.064 - 0.128' },
  'J': { altitude: '28587 - 190580', weight: '3.2 - 6.4', payload: '0.128 - 0.256' },
  'K': { altitude: '114348 - 762320', weight: '6.4 - 12.8', payload: '0.256 - 0.512' },
  'L': { altitude: '457863 - 3.05e6', weight: '12.8 - 25.6', payload: '0.512 - 1.024' },
  'M': { altitude: '1.83e6 - 1.22e7', weight: '25.6 - 51.2', payload: '1.024 - 2.048' },
  'N': { altitude: '7.32e6 - 4.88e7', weight: '51.2 - 102.4', payload: '2.048 - 4.096' },
  'O': { altitude: '2.93e7 - 1.95e8', weight: '102.4 - 204.8', payload: '4.096 - 8.192' },
  'P': { altitude: '1.17e8 - 7.79e8', weight: '204.8 - 409.6', payload: '8.192 - 16.384' },
  'Q': { altitude: '4.68e8 - 3.12e9', weight: '409.6 - 819.2', payload: '16.384 - 32.768' },
  'R': { altitude: '1.87e9 - 1.25e10', weight: '819.2 - 1638.4', payload: '32.768 - 65.536' },
  'S': { altitude: '7.49e9 - 4.99e10', weight: '1638.4 - 3276.8', payload: '65.536 - 131.072' },
  'T': { altitude: '2.99e10 - 1.99e11', weight: '3276.8 - 6553.6', payload: '131.072 - 262.144' },
  'U': { altitude: '1.20e11 - 7.97e11', weight: '6553.6 - 13107.2', payload: '262.144 - 524.288' },
  'V': { altitude: '4.79e11 - 3.19e12', weight: '13107.2 - 26214.4', payload: '524.288 - 1048.576' }
};

const App: React.FC = () => {
  // Authentication State
  const [isAuthenticated, setIsAuthenticated] = useState(false);
  const [username, setUsername] = useState('zulfam');
  const [password, setPassword] = useState('');
  const [authError, setAuthError] = useState('');

  // App State
  const [inputs, setInputs] = useState<MotorInputs>(INITIAL_INPUTS);
  const [results, setResults] = useState<SimulationResults | null>(null);
  const [activeTab, setActiveTab] = useState<'data' | 'plots'>('data');

  const handleLogin = (e: React.FormEvent) => {
    e.preventDefault();
    if (username === 'zulfamadnan' && password === 'Zpassword@123') {
      setIsAuthenticated(true);
      setAuthError('');
    } else {
      setAuthError('Invalid username or password');
    }
  };

  const handleLogout = () => {
    setIsAuthenticated(false);
    setPassword('');
    setAuthError('');
    setResults(null);
  };

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
    const { name, value, type } = e.target;
    setInputs(prev => ({
      ...prev,
      [name]: type === 'number' ? parseFloat(value) : value
    }));
  };

  const currentExitDiameter = useMemo(() => {
     return inputs.throatDiameter * Math.sqrt(inputs.nozzleExpansionRatio);
  }, [inputs.throatDiameter, inputs.nozzleExpansionRatio]);

  const calculateSimulation = () => {
    // 1. Calculate Constants
    const PI = Math.PI;
    const { 
      propellantOD, 
      propellantCore, 
      propellantLength, 
      chamberDiameter, 
      totalLength, 
      throatDiameter,
      throatErosion,
      exposureOption,
      burnRateA,
      burnRateN,
      molecularMass,
      temperature,
      specificHeatRatio,
      density,
      atmosphericPressure,
      nozzleExpansionRatio,
      nozzleEfficiency,
      propellantGStar
    } = inputs;

    // --- Thermodynamics Calculations ---
    // R = 8314 / MolecularMass
    const R = 8314 / molecularMass;

    // C* Calculation
    const k = specificHeatRatio;
    const term1 = (R * temperature) / k;
    const term2Base = (k + 1) / 2;
    const term2Exp = (k + 1) / (k - 1);
    const term2 = Math.pow(term2Base, term2Exp);
    const cStar = Math.sqrt(term1 * term2);

    // Convert dimensions to consistent float values
    const D_out = propellantOD;
    const D_core = propellantCore;
    const L_grain = propellantLength;
    const D_throat_initial = throatDiameter;

    // --- Nozzle Geometry Calculations ---
    const nozzleExitDiameter = D_throat_initial * Math.sqrt(nozzleExpansionRatio);

    // Basic Volumes
    const r_chamber = chamberDiameter / 2;
    const chamberVol = PI * Math.pow(r_chamber, 2) * totalLength;

    const r_out = D_out / 2;
    const r_core = D_core / 2;
    const propCrossSectionArea = PI * (Math.pow(r_out, 2) - Math.pow(r_core, 2));
    const propellantVol = propCrossSectionArea * L_grain;
    
    // Mass = Density * Volume
    // Input Density is in g/cm³ -> kg/m³ conversion
    const mm3_to_m3 = 1e-9;
    const mm2_to_m2 = 1e-6;
    const g_cm3_to_kg_m3 = 1000;
    const rho_prop_kgm3 = density * g_cm3_to_kg_m3;
    
    // Mass (g) = (Density (g/cm³) / 1000) * Volume (mm³)
    const propellantMass = (density / 1000) * propellantVol;

    // Web Thickness (Radius difference)
    const webThickness = (D_out - D_core) / 2;

    // Simulation Loop Setup
    const steps = 1000;
    const stepSize = webThickness / steps; // mm
    const data: SimulationStep[] = [];

    // Physics Simulation State
    let t = 0; // seconds
    let P_abs = atmosphericPressure; // MPa (Initial = Patm)
    let M_stored = 0; // kg
    let P_max = P_abs;
    let F_max = 0;
    let TotalImpulseBurn = 0;

    const chamberVol_m3 = chamberVol * mm3_to_m3;
    const P_atm_Pa = atmosphericPressure * 1e6;

    // Helper: Thrust Calculation
    const calculateThrust = (Pi_MPa: number, At_mm2: number) => {
        if (Pi_MPa <= atmosphericPressure) return 0;

        const Pi_Pa = Pi_MPa * 1e6;
        const At_m2 = At_mm2 * mm2_to_m2;
        
        const k2 = k * k;
        const exponent1 = (k + 1) / (k - 1);
        const exponent2 = (k - 1) / k;
        
        const pressureRatio = P_atm_Pa / Pi_Pa;
        if (pressureRatio >= 1) return 0;

        const bracketTerm = 1 - Math.pow(pressureRatio, exponent2);
        const coeffTerm = (2 * k2) / (k - 1) * Math.pow(2 / (k + 1), exponent1);
        
        const Cf = Math.sqrt(coeffTerm * bracketTerm);
        const F = Pi_Pa * At_m2 * nozzleEfficiency * Cf;
        return F > 0 ? F : 0;
    };

    // Initial Geometry State
    let currentID = D_core;
    let currentThroatDiameter = D_throat_initial;
    
    // Calculate initial V_prop for mass diff
    const getPropVol = (id: number) => {
        const r_o = D_out / 2;
        const r_i = id / 2;
        if (id >= D_out) return 0;
        const area = PI * (Math.pow(r_o, 2) - Math.pow(r_i, 2));
        return area * L_grain; // mm3
    };

    let initialKn = 0;
    let maxKn = 0;
    let initialAb = 0;
    let finalThroatDiameter = D_throat_initial;
    let stepCounter = 0;

    // Step 0 Data
    {
        const Area_port = PI * currentID * L_grain;
        let Area_ends = 0;
        if (exposureOption === ExposureOption.PORT_AND_ENDS) {
             Area_ends = 2 * (PI / 4) * (Math.pow(D_out, 2) - Math.pow(currentID, 2));
        }
        const totalBurningArea = Area_port + Area_ends;
        const currentA_throat = PI * Math.pow(currentThroatDiameter / 2, 2);
        const kn = currentA_throat > 0 ? totalBurningArea / currentA_throat : 0;
        
        initialKn = kn;
        initialAb = totalBurningArea;
        maxKn = kn;

        const thrust = calculateThrust(P_abs, currentA_throat);

        data.push({
            step: stepCounter++,
            webRegression: 0,
            webRemaining: webThickness,
            kn: parseFloat(kn.toFixed(2)),
            burningArea: parseFloat(totalBurningArea.toFixed(1)),
            throatDiameter: parseFloat(currentThroatDiameter.toFixed(2)),
            pressure: parseFloat(P_abs.toFixed(4)),
            thrust: parseFloat(thrust.toFixed(2)),
            time: 0
        });
    }

    let last_V_free_m3 = 0;

    // --- MAIN INTEGRATION LOOP ---
    for (let i = 1; i <= steps; i++) {
      const x = i * stepSize;
      const x_prev = (i - 1) * stepSize;
      
      const currID = D_core + (2 * x);
      const prevID = D_core + (2 * x_prev);

      if (currID > D_out + 0.001) break;

      const P_prev = P_abs;

      // 1. Calculate mn (Mass Flow Nozzle) based on P(i-1)
      // Uses throat area at start of step
      const erosionProgPrev = webThickness > 0 ? (x_prev / webThickness) : 1;
      const Dt_prev = D_throat_initial + (throatErosion * erosionProgPrev);
      const At_prev_m2 = PI * Math.pow(Dt_prev / 2, 2) * mm2_to_m2;

      const P_gauge_Pa_prev = (P_prev - atmosphericPressure) * 1e6;
      let m_n = 0; // kg/s

      if (P_gauge_Pa_prev > 0) {
         // mn formula: ((Pi-1)-patm)*1e6 * At / SQRT(R*T) * SQRT(k) * (2/(k+1))^((k+1)/2/(k-1))
         const num = P_gauge_Pa_prev * At_prev_m2;
         const den = Math.sqrt(R * temperature);
         const k_term = Math.sqrt(k) * Math.pow(2 / (k + 1), (k + 1) / (2 * (k - 1)));
         m_n = (num / den) * k_term;
      }

      // 2. Geometry for Erosive Calculation (at start of step)
      const Dp_m = prevID / 1000; // Dpi
      
      // Calculate current propellant length (Lpi)
      let currentLen_mm = propellantLength;
      if (exposureOption === ExposureOption.PORT_AND_ENDS) {
          // If ends are burning, length reduces by 2*x
          currentLen_mm = propellantLength - (2 * x_prev);
      }
      if (currentLen_mm < 1) currentLen_mm = 1;
      const Lp_m = currentLen_mm / 1000; // Lpi in meters

      const Ap_m2 = PI * Math.pow(Dp_m / 2, 2); // Api
      const pp_m = PI * Dp_m; // ppi

      // 3. Mass Flux M (Mi)
      const M = Ap_m2 > 0 ? m_n / Ap_m2 : 0;

      // 4. Burn Rate Calculation (rei)
      // rei = (1 + ((Lpi/Dpi)/(3.14*T )) * ((P(i-1))^-0.7) * ((3.14*Api/ppi)^-0.2) * (Mi^0.65*(1-(10/Mi)))) * aPi^n
      const r_base = burnRateA * Math.pow(P_prev, burnRateN);
      let r = r_base;

      const threshold = 10; // Explicitly 10 as per new formula

      if (M > threshold && P_prev > 1e-4) {
         // Term 1: (Lpi/Dpi)/(3.14*T)
         const term1 = (Lp_m / Dp_m) / (3.14 * temperature);
         
         // Term 2: (P(i-1))^-0.7
         const term2 = Math.pow(P_prev, -0.7);
         
         // Term 3: ((3.14*Api/ppi)^-0.2)
         const term3_inner = 3.14 * (Ap_m2 / pp_m);
         const term3 = Math.pow(term3_inner, -0.2);
         
         // Term 4: (Mi^0.65*(1-(10/Mi)))
         const term4 = Math.pow(M, 0.65) * (1 - (threshold / M));

         const erosiveFactor = term1 * term2 * term3 * term4;
         
         r = r_base * (1 + erosiveFactor);
      }

      // 5. Time Step Evaluation
      const safe_r = r > 0.001 ? r : 0.001; 
      const dt = stepSize / safe_r;
      t += dt; // Update time evaluated at each step

      // 6. Mass Generation (in Chamber)
      const V_prop_prev_mm3 = getPropVol(prevID);
      const V_prop_curr_mm3 = getPropVol(currID);
      const V_burned_mm3 = V_prop_prev_mm3 - V_prop_curr_mm3;
      const dm_gen = (V_burned_mm3 * mm3_to_m3) * rho_prop_kgm3;

      // 7. Mass Out (from Nozzle) over dt
      const dm_out = m_n * dt;

      // 8. Update Mass Stored
      const dm_stored = dm_gen - dm_out;
      M_stored += dm_stored; // Update mass stored in chamber
      if (M_stored < 0) M_stored = 0;

      // 9. Evaluate Chamber Pressure (P_abs) at each step
      const V_free_m3 = chamberVol_m3 - (V_prop_curr_mm3 * mm3_to_m3);
      last_V_free_m3 = V_free_m3;

      const rho_stored = V_free_m3 > 0 ? M_stored / V_free_m3 : 0;
      const P_calc_gauge_Pa = rho_stored * R * temperature;
      P_abs = (P_calc_gauge_Pa / 1e6) + atmosphericPressure;

      if (P_abs > P_max) P_max = P_abs;

      // --- END OF INTEGRATION STEP ---

      // Calculate auxiliary outputs
      const erosionProgressCurrent = webThickness > 0 ? (x / webThickness) : 1;
      const currentThroatDia = D_throat_initial + (throatErosion * erosionProgressCurrent);
      finalThroatDiameter = currentThroatDia;

      const currentA_throat_curr = PI * Math.pow(currentThroatDia / 2, 2);
      
      const Area_port = PI * currID * L_grain;
      let Area_ends = 0;
      if (exposureOption === ExposureOption.PORT_AND_ENDS) {
        Area_ends = 2 * (PI / 4) * (Math.pow(D_out, 2) - Math.pow(currID, 2));
      }
      const totalBurningArea = Area_port + Area_ends;
      const kn = currentA_throat_curr > 0 ? totalBurningArea / currentA_throat_curr : 0;
      if (kn > maxKn) maxKn = kn;

      const F_curr = calculateThrust(P_abs, currentA_throat_curr);
      if (F_curr > F_max) F_max = F_curr;
      
      TotalImpulseBurn += F_curr * dt;

      data.push({
        step: stepCounter++,
        webRegression: parseFloat(x.toFixed(3)),
        webRemaining: parseFloat((webThickness - x).toFixed(3)),
        kn: parseFloat(kn.toFixed(2)),
        burningArea: parseFloat(totalBurningArea.toFixed(1)),
        throatDiameter: parseFloat(currentThroatDia.toFixed(2)),
        pressure: parseFloat(P_abs.toFixed(4)),
        thrust: parseFloat(F_curr.toFixed(2)),
        time: parseFloat(t.toFixed(4))
      });
    }

    const calculatedBurnOutTime = t;
    const F_avg = calculatedBurnOutTime > 0 ? TotalImpulseBurn / calculatedBurnOutTime : 0;
    
    // Motor Class Calculation
    const motorClass = getMotorClass(TotalImpulseBurn);
    
    // Estimates based on class
    const estimates = CLASS_ESTIMATES[motorClass] || { altitude: '-', weight: '-', payload: '-' };

    // --- Specific Impulse Calculation ---
    const propellantMassKg = propellantMass / 1000;
    const gravity = 9.8;
    let specificImpulse = 0;

    if (calculatedBurnOutTime > 0 && propellantMassKg > 0) {
        const massFlowRate = propellantMassKg / calculatedBurnOutTime;
        specificImpulse = F_avg / (massFlowRate * gravity);
    }

    // --- Tail-off Simulation ---
    const tailOffDt = 0.005;
    const maxTailOffSteps = 500;
    let tailOffStep = 0;
    const V_free_tailoff_m3 = last_V_free_m3 > 0 ? last_V_free_m3 : chamberVol_m3;

    while (P_abs > atmosphericPressure * 1.01 && tailOffStep < maxTailOffSteps) {
        tailOffStep++;
        t += tailOffDt;

        const P_gauge_MPa = P_abs - atmosphericPressure;
        const P_drive_Pa = P_gauge_MPa * 1e6;
        
        // Flow calc
        const A_throat_curr = PI * Math.pow(finalThroatDiameter / 2, 2);
        const A_throat_m2 = A_throat_curr * mm2_to_m2;
        let m_out_rate = 0;

        if (P_drive_Pa > 0) {
            const denomTerm1 = Math.sqrt(R * temperature);
            const k_term = Math.sqrt(k) * Math.pow(2 / (k + 1), (k + 1) / (2 * (k - 1)));
            m_out_rate = (P_drive_Pa * A_throat_m2 / denomTerm1) * k_term;
        }

        const dm_stored = (0 - m_out_rate) * tailOffDt;
        M_stored += dm_stored;
        if (M_stored < 0) M_stored = 0;

        const rho_stored = M_stored / V_free_tailoff_m3;
        const P_calc_gauge_Pa = rho_stored * R * temperature;
        P_abs = (P_calc_gauge_Pa / 1e6) + atmosphericPressure;
        
        // Thrust Calc
        const F_curr = calculateThrust(P_abs, A_throat_curr);

        data.push({
            step: stepCounter++,
            webRegression: parseFloat(webThickness.toFixed(3)),
            webRemaining: 0,
            kn: 0,
            burningArea: 0,
            throatDiameter: parseFloat(finalThroatDiameter.toFixed(2)),
            pressure: parseFloat(P_abs.toFixed(4)),
            thrust: parseFloat(F_curr.toFixed(2)),
            time: parseFloat(t.toFixed(4))
        });
    }

    // Burn Rate Curve
    const burnRateData: BurnRateStep[] = [];
    for (let p = 0; p <= 10; p += 0.2) {
      const r = burnRateA * Math.pow(p, burnRateN);
      burnRateData.push({
        pressure: parseFloat(p.toFixed(1)),
        rate: parseFloat(r.toFixed(2))
      });
    }

    setResults({
      chamberVolume: parseFloat(chamberVol.toFixed(2)),
      propellantVolume: parseFloat(propellantVol.toFixed(2)),
      propellantMass: parseFloat(propellantMass.toFixed(1)),
      initialKn: parseFloat(initialKn.toFixed(2)),
      maxKn: parseFloat(maxKn.toFixed(2)),
      maxPressure: parseFloat(P_max.toFixed(2)),
      maxThrust: parseFloat(F_max.toFixed(1)),
      averageThrust: parseFloat(F_avg.toFixed(1)),
      totalImpulse: parseFloat(TotalImpulseBurn.toFixed(1)),
      motorClass: motorClass,
      estimatedAltitude: estimates.altitude,
      estimatedRocketWeight: estimates.weight,
      estimatedPayload: estimates.payload,
      specificImpulse: parseFloat(specificImpulse.toFixed(1)),
      burnTime: parseFloat(t.toFixed(2)),
      burnOutTime: parseFloat(calculatedBurnOutTime.toFixed(2)),
      webThickness: parseFloat(webThickness.toFixed(2)),
      initialBurningArea: parseFloat(initialAb.toFixed(2)),
      finalThroatDiameter: parseFloat(finalThroatDiameter.toFixed(2)),
      nozzleExitDiameter: parseFloat(nozzleExitDiameter.toFixed(2)),
      gasConstantR: parseFloat(R.toFixed(1)),
      cStar: parseFloat(cStar.toFixed(1)),
      data,
      burnRateData
    });
  };

  if (!isAuthenticated) {
    return (
      <div className="min-h-screen bg-slate-50 flex items-center justify-center p-4">
        <div className="bg-white border border-slate-200 p-8 rounded-xl shadow-xl w-full max-w-md">
          <div className="flex justify-center mb-6">
            <div className="bg-indigo-600 p-3 rounded-lg shadow-lg shadow-indigo-500/20">
              <Rocket className="w-8 h-8 text-white" />
            </div>
          </div>
          <h2 className="text-2xl font-bold text-center text-slate-900 mb-2">Simulator Login</h2>
          <p className="text-slate-500 text-center text-sm mb-8">Enter your credentials to access the tool</p>
          
          <form onSubmit={handleLogin} className="space-y-4">
            <div>
              <label className="block text-xs font-bold text-slate-600 uppercase tracking-wider mb-1.5 ml-1">Username</label>
              <div className="relative">
                <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                  <User className="h-4 w-4 text-slate-400" />
                </div>
                <input 
                  type="text" 
                  value={username}
                  onChange={(e) => setUsername(e.target.value)}
                  className="w-full bg-slate-50 border border-slate-300 rounded-lg py-2.5 pl-10 pr-3 text-slate-900 placeholder-slate-400 focus:ring-2 focus:ring-indigo-500 focus:border-transparent outline-none transition-all"
                  placeholder="Enter username"
                />
              </div>
            </div>
            
            <div>
              <label className="block text-xs font-bold text-slate-600 uppercase tracking-wider mb-1.5 ml-1">Password</label>
              <div className="relative">
                <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                  <Lock className="h-4 w-4 text-slate-400" />
                </div>
                <input 
                  type="password" 
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                  className="w-full bg-slate-50 border border-slate-300 rounded-lg py-2.5 pl-10 pr-3 text-slate-900 placeholder-slate-400 focus:ring-2 focus:ring-indigo-500 focus:border-transparent outline-none transition-all"
                  placeholder="Enter password"
                />
              </div>
            </div>
            
            {authError && (
              <div className="bg-red-50 border border-red-200 text-red-600 text-xs p-3 rounded-lg text-center font-medium">
                {authError}
              </div>
            )}

            <button 
              type="submit"
              className="w-full bg-indigo-600 hover:bg-indigo-700 text-white font-bold py-3 px-6 rounded-lg shadow-lg shadow-indigo-600/20 transition-all hover:scale-[1.02] active:scale-[0.98] mt-2 flex items-center justify-center space-x-2"
            >
              <span>Sign In</span>
            </button>
          </form>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900 p-4 md:p-8 font-sans">
      <header className="max-w-7xl mx-auto mb-8 flex flex-col md:flex-row md:items-center justify-between border-b border-slate-200 pb-6 gap-4">
        <div className="flex items-center space-x-3">
          <div className="bg-indigo-600 p-2.5 rounded-lg shadow-lg shadow-indigo-500/20">
            <Rocket className="w-8 h-8 text-white" />
          </div>
          <div>
            <h1 className="text-3xl font-bold tracking-tight text-slate-900">SRM Simulator</h1>
            <p className="text-slate-500 text-sm">Solid Rocket Motor Ballistic Analysis Tool</p>
          </div>
        </div>
        
        <div className="flex items-center gap-4">
          <div className="text-right hidden md:block">
            <p className="text-sm font-medium text-slate-900">{username}</p>
            <p className="text-xs text-slate-500">Authorized User</p>
          </div>
          <button 
            onClick={handleLogout}
            className="flex items-center space-x-2 px-4 py-2 bg-white hover:bg-slate-50 border border-slate-200 rounded-lg text-slate-600 hover:text-slate-900 transition-all text-sm font-medium shadow-sm"
          >
            <LogOut className="w-4 h-4" />
            <span>Logout</span>
          </button>
        </div>
      </header>

      {/* TABS NAVIGATION */}
      <div className="max-w-7xl mx-auto mb-8">
        <div className="flex space-x-2 bg-white p-1 rounded-xl border border-slate-200 w-fit shadow-sm">
          <button
            onClick={() => setActiveTab('data')}
            className={`flex items-center space-x-2 px-6 py-2.5 rounded-lg text-sm font-medium transition-all ${
              activeTab === 'data' 
                ? 'bg-indigo-600 text-white shadow-md shadow-indigo-500/30' 
                : 'text-slate-600 hover:text-slate-900 hover:bg-slate-100'
            }`}
          >
            <LayoutDashboard className="w-4 h-4" />
            <span>Data & Configuration</span>
          </button>
          <button
            onClick={() => setActiveTab('plots')}
            className={`flex items-center space-x-2 px-6 py-2.5 rounded-lg text-sm font-medium transition-all ${
              activeTab === 'plots' 
                ? 'bg-indigo-600 text-white shadow-md shadow-indigo-500/30' 
                : 'text-slate-600 hover:text-slate-900 hover:bg-slate-100'
            }`}
          >
            <LineChartIcon className="w-4 h-4" />
            <span>Simulation Plots</span>
          </button>
        </div>
      </div>

      <main className="max-w-7xl mx-auto">
        {activeTab === 'data' ? (
          <div className="grid grid-cols-1 lg:grid-cols-12 gap-8">
            {/* INPUT PANEL - White Background */}
            <section className="lg:col-span-4 space-y-6">
              <div className="bg-white border border-slate-200 rounded-xl p-6 shadow-sm">
                <div className="flex items-center space-x-2 mb-6 text-indigo-600">
                  <Ruler className="w-5 h-5" />
                  <h2 className="text-lg font-bold text-slate-900">Motor Geometry</h2>
                </div>
                <div className="grid grid-cols-2 gap-4">
                  <InputGroup label="Total Length" name="totalLength" unit="mm" value={inputs.totalLength} onChange={handleInputChange} />
                  <InputGroup label="Head/Nose Len" name="headLength" unit="mm" value={inputs.headLength} onChange={handleInputChange} />
                  <InputGroup label="Chamber ID" name="chamberDiameter" unit="mm" value={inputs.chamberDiameter} onChange={handleInputChange} />
                  <InputGroup label="Propellant OD" name="propellantOD" unit="mm" value={inputs.propellantOD} onChange={handleInputChange} />
                  <InputGroup label="Core Diameter" name="propellantCore" unit="mm" value={inputs.propellantCore} onChange={handleInputChange} />
                  <div className="col-span-2">
                    <InputGroup label="Propellant Length" name="propellantLength" unit="mm" value={inputs.propellantLength} onChange={handleInputChange} />
                  </div>
                </div>
              </div>

              <div className="bg-white border border-slate-200 rounded-xl p-6 shadow-sm">
                <div className="flex items-center space-x-2 mb-6 text-orange-500">
                  <Flame className="w-5 h-5" />
                  <h2 className="text-lg font-bold text-slate-900">Propellant & Nozzle</h2>
                </div>
                <div className="space-y-4">
                   <InputGroup 
                      label="Burning Surfaces" 
                      name="exposureOption" 
                      type="select" 
                      options={[ExposureOption.PORT_ONLY, ExposureOption.PORT_AND_ENDS]} 
                      value={inputs.exposureOption} 
                      onChange={handleInputChange} 
                    />
                  <div className="grid grid-cols-2 gap-4">
                     <InputGroup label="Density" name="density" unit="g/cm³" step="0.01" value={inputs.density} onChange={handleInputChange} />
                     <InputGroup label="Throat Ø" name="throatDiameter" unit="mm" value={inputs.throatDiameter} onChange={handleInputChange} />
                  </div>
                  
                  <div className="border-t border-slate-100 pt-4 mt-2">
                     <div className="flex items-center space-x-1 mb-3">
                       <Cone className="w-4 h-4 text-slate-400" />
                       <h3 className="text-xs font-bold text-slate-500 uppercase tracking-wider">Nozzle Design</h3>
                     </div>
                     <div className="grid grid-cols-2 gap-4">
                        <InputGroup label="Conv. Angle (β)" name="nozzleConvergingAngle" unit="deg" value={inputs.nozzleConvergingAngle} onChange={handleInputChange} />
                        <InputGroup label="Div. Angle (α)" name="nozzleDivergingAngle" unit="deg" value={inputs.nozzleDivergingAngle} onChange={handleInputChange} />
                        <InputGroup label="Expansion Ratio" name="nozzleExpansionRatio" unit="Ae/At" step="0.1" value={inputs.nozzleExpansionRatio} onChange={handleInputChange} />
                        <InputGroup label="Efficiency (Nef)" name="nozzleEfficiency" unit="0-1" step="0.01" value={inputs.nozzleEfficiency} onChange={handleInputChange} />
                     </div>
                  </div>

                  <div>
                     <InputGroup 
                       label="Total Throat Erosion" 
                       name="throatErosion" 
                       unit="mm" 
                       description="Total diameter increase by burnout"
                       value={inputs.throatErosion} 
                       onChange={handleInputChange} 
                     />
                  </div>
                </div>
              </div>

              <div className="bg-white border border-slate-200 rounded-xl p-6 shadow-sm">
                 <div className="flex items-center space-x-2 mb-6 text-emerald-600">
                  <Settings className="w-5 h-5" />
                  <h2 className="text-lg font-bold text-slate-900">Thermodynamics</h2>
                </div>
                <div className="grid grid-cols-2 gap-4">
                   <InputGroup label="Specific Heat (k)" name="specificHeatRatio" unit="const" step="0.01" value={inputs.specificHeatRatio} onChange={handleInputChange} />
                   <InputGroup label="Temperature" name="temperature" unit="K" value={inputs.temperature} onChange={handleInputChange} />
                   <InputGroup label="Mol. Mass" name="molecularMass" unit="g/mol" value={inputs.molecularMass} onChange={handleInputChange} />
                   <InputGroup label="G* Threshold" name="propellantGStar" unit="g/s-cm²" value={inputs.propellantGStar} onChange={handleInputChange} />
                   <InputGroup label="Atm Pressure" name="atmosphericPressure" unit="MPa" step="0.000001" value={inputs.atmosphericPressure} onChange={handleInputChange} />
                   <InputGroup label="Coeff (a)" name="burnRateA" unit="const" value={inputs.burnRateA} onChange={handleInputChange} />
                   <div className="col-span-2">
                     <InputGroup label="Exponent (n)" name="burnRateN" unit="const" value={inputs.burnRateN} onChange={handleInputChange} />
                   </div>
                </div>
              </div>

              <button 
                onClick={calculateSimulation}
                className="w-full bg-indigo-600 hover:bg-indigo-700 text-white font-bold py-4 px-6 rounded-xl flex items-center justify-center space-x-2 transition-all transform active:scale-95 shadow-lg shadow-indigo-600/30"
              >
                <Play className="w-5 h-5 fill-current" />
                <span>Evaluate Motor</span>
              </button>
            </section>

            {/* OUTPUT PANEL - Distinct Gray Background or Cards */}
            <section className="lg:col-span-8 space-y-6">
              {/* Visualizer */}
              <MotorVisualizer inputs={inputs} calculatedExitDiameter={currentExitDiameter} />

              {/* Output Cards Row 1: Geometric Props */}
              <div className="grid grid-cols-1 md:grid-cols-3 xl:grid-cols-5 gap-4">
                 <ResultCard 
                    label="Chamber Vol" 
                    value={results ? (results.chamberVolume / 1000).toFixed(1) : '-'} 
                    unit="cc" 
                    icon={Cylinder} 
                 />
                 <ResultCard 
                    label="Propellant Vol" 
                    value={results ? (results.propellantVolume / 1000).toFixed(1) : '-'} 
                    unit="cc"
                    icon={Box} 
                 />
                 <ResultCard 
                    label="Propellant Mass" 
                    value={results ? results.propellantMass : '-'} 
                    unit="g" 
                    icon={Weight} 
                 />
                 <ResultCard 
                    label="Initial Area" 
                    value={results ? (results.initialBurningArea / 100).toFixed(1) : '-'} 
                    unit="cm²"
                    icon={Maximize} 
                 />
                 <ResultCard 
                    label="Nozzle Exit Ø" 
                    value={results ? results.nozzleExitDiameter : '-'} 
                    unit="mm" 
                    icon={Cone} 
                 />
              </div>

              {/* Output Cards Row 2: Performance Metrics */}
              <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-4 gap-4">
                 <ResultCard 
                    label="Motor Class" 
                    value={results ? results.motorClass : '-'} 
                    unit=""
                    icon={Award} 
                 />
                 <ResultCard 
                    label="Total Impulse" 
                    value={results ? results.totalImpulse : '-'} 
                    unit="Ns"
                    icon={Zap} 
                 />
                 <ResultCard 
                    label="Max Thrust" 
                    value={results ? results.maxThrust : '-'} 
                    unit="N"
                    icon={ArrowUp} 
                 />
                 <ResultCard 
                    label="Avg Thrust" 
                    value={results ? results.averageThrust : '-'} 
                    unit="N"
                    subValue={results ? "During Burn" : undefined}
                    icon={TrendingUp} 
                 />
                 {/* New Estimates */}
                 <ResultCard 
                    label="Est. Altitude" 
                    value={results ? results.estimatedAltitude : '-'} 
                    unit="m"
                    icon={Mountain} 
                 />
                 <ResultCard 
                    label="Est. Rocket Wt" 
                    value={results ? results.estimatedRocketWeight : '-'} 
                    unit="kg"
                    icon={Scale} 
                 />
                 <ResultCard 
                    label="Est. Payload" 
                    value={results ? results.estimatedPayload : '-'} 
                    unit="kg"
                    icon={Package} 
                 />
                 <ResultCard 
                    label="Specific Impulse" 
                    value={results ? results.specificImpulse : '-'} 
                    unit="s"
                    icon={Activity} 
                 />
                 <ResultCard 
                    label="Max Pressure" 
                    value={results ? results.maxPressure : '-'} 
                    unit="MPa"
                    icon={Gauge} 
                 />
                 <ResultCard 
                    label="Burn Out Time" 
                    value={results ? results.burnOutTime : '-'} 
                    unit="s"
                    icon={Flame} 
                 />
                 <ResultCard 
                    label="Total Time" 
                    value={results ? results.burnTime : '-'} 
                    unit="s"
                    icon={Timer} 
                 />
              </div>
            </section>
          </div>
        ) : (
          /* PLOTS TAB */
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
             {/* Chart 1: Time vs Thrust - Output Panel Style */}
             <div className="bg-slate-100 border border-slate-200 rounded-xl p-6 shadow-sm h-[450px] flex flex-col lg:col-span-2">
              <div className="flex justify-between items-center mb-6">
                <div className="flex items-center space-x-2 text-cyan-600">
                  <ArrowUp className="w-5 h-5" />
                  <h3 className="text-lg font-bold text-slate-900">Thrust Curve</h3>
                </div>
                {results && (
                  <div className="flex space-x-4 text-xs font-mono text-slate-500">
                     <span>Fmax: <span className="text-slate-900 font-bold">{results.maxThrust} N</span></span>
                     <span>Isp: <span className="text-slate-900 font-bold">{results.specificImpulse} s</span></span>
                  </div>
                )}
              </div>

              <div className="flex-1 w-full min-h-0">
                  {results ? (
                    <ResponsiveContainer width="100%" height="100%">
                      <AreaChart data={results.data} margin={{ top: 10, right: 30, left: 10, bottom: 20 }}>
                        <defs>
                          <linearGradient id="colorThrust" x1="0" y1="0" x2="0" y2="1">
                            <stop offset="5%" stopColor="#0891b2" stopOpacity={0.4}/>
                            <stop offset="95%" stopColor="#0891b2" stopOpacity={0}/>
                          </linearGradient>
                        </defs>
                        <CartesianGrid strokeDasharray="3 3" stroke="#cbd5e1" opacity={0.6} />
                        <XAxis 
                          dataKey="time" 
                          stroke="#64748b" 
                          label={{ value: 'Time (s)', position: 'insideBottom', offset: -10, fill: '#64748b' }} 
                          tick={{fontSize: 12}}
                          allowDecimals={true}
                          tickFormatter={(val) => val.toFixed(2)}
                        />
                        <YAxis 
                          stroke="#0891b2" 
                          label={{ value: 'Thrust (N)', angle: -90, position: 'insideLeft', fill: '#0891b2', offset: 0 }}
                          tick={{fontSize: 12}}
                        />
                        <Tooltip 
                          contentStyle={{ backgroundColor: '#ffffff', borderColor: '#e2e8f0', borderRadius: '8px', color: '#0f172a' }}
                          formatter={(value: number) => [value, 'N']}
                          labelFormatter={(label) => `Time: ${parseFloat(label).toFixed(3)} s`}
                        />
                        <Area 
                          type="monotone" 
                          dataKey="thrust" 
                          name="Thrust"
                          stroke="#0891b2" 
                          strokeWidth={3} 
                          fillOpacity={1} 
                          fill="url(#colorThrust)" 
                          animationDuration={1500}
                        />
                      </AreaChart>
                    </ResponsiveContainer>
                  ) : (
                    <div className="h-full flex items-center justify-center text-slate-500">Run evaluation to view data</div>
                  )}
              </div>
            </div>

            {/* Chart 2: Time vs Pressure */}
            <div className="bg-slate-100 border border-slate-200 rounded-xl p-6 shadow-sm h-[450px] flex flex-col">
              <div className="flex justify-between items-center mb-6">
                <div className="flex items-center space-x-2 text-rose-500">
                  <TrendingUp className="w-5 h-5" />
                  <h3 className="text-lg font-bold text-slate-900">Chamber Pressure</h3>
                </div>
                {results && (
                  <div className="flex space-x-4 text-xs font-mono text-slate-500">
                     <span>Max P: <span className="text-slate-900 font-bold">{results.maxPressure} MPa</span></span>
                  </div>
                )}
              </div>

              <div className="flex-1 w-full min-h-0">
                  {results ? (
                    <ResponsiveContainer width="100%" height="100%">
                      <AreaChart data={results.data} margin={{ top: 10, right: 30, left: 10, bottom: 20 }}>
                        <defs>
                          <linearGradient id="colorPressure" x1="0" y1="0" x2="0" y2="1">
                            <stop offset="5%" stopColor="#e11d48" stopOpacity={0.4}/>
                            <stop offset="95%" stopColor="#e11d48" stopOpacity={0}/>
                          </linearGradient>
                        </defs>
                        <CartesianGrid strokeDasharray="3 3" stroke="#cbd5e1" opacity={0.6} />
                        <XAxis 
                          dataKey="time" 
                          stroke="#64748b" 
                          label={{ value: 'Time (s)', position: 'insideBottom', offset: -10, fill: '#64748b' }} 
                          tick={{fontSize: 12}}
                          allowDecimals={true}
                          tickFormatter={(val) => val.toFixed(2)}
                        />
                        <YAxis 
                          stroke="#e11d48" 
                          label={{ value: 'Pressure (MPa)', angle: -90, position: 'insideLeft', fill: '#e11d48', offset: 0 }}
                          tick={{fontSize: 12}}
                        />
                        <Tooltip 
                          contentStyle={{ backgroundColor: '#ffffff', borderColor: '#e2e8f0', borderRadius: '8px', color: '#0f172a' }}
                          formatter={(value: number) => [value, 'MPa']}
                          labelFormatter={(label) => `Time: ${parseFloat(label).toFixed(3)} s`}
                        />
                        <Area 
                          type="monotone" 
                          dataKey="pressure" 
                          name="Chamber Pressure"
                          stroke="#e11d48" 
                          strokeWidth={3} 
                          fillOpacity={1} 
                          fill="url(#colorPressure)" 
                          animationDuration={1500}
                        />
                      </AreaChart>
                    </ResponsiveContainer>
                  ) : (
                    <div className="h-full flex items-center justify-center text-slate-500">Run evaluation to view data</div>
                  )}
              </div>
            </div>

            {/* Chart 3: Kn vs Web Regression */}
            <div className="bg-slate-100 border border-slate-200 rounded-xl p-6 shadow-sm h-[450px] flex flex-col">
              <div className="flex justify-between items-center mb-6">
                <div className="flex items-center space-x-2 text-indigo-600">
                  <Activity className="w-5 h-5" />
                  <h3 className="text-lg font-bold text-slate-900">Web Regression</h3>
                </div>
                {results && (
                   <div className="flex space-x-4 text-xs font-mono text-slate-500">
                      <span>Max Kn: <span className="text-slate-900 font-bold">{results.maxKn}</span></span>
                   </div>
                )}
              </div>

              <div className="flex-1 w-full min-h-0">
                 {results ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <ComposedChart data={results.data} margin={{ top: 10, right: 30, left: 10, bottom: 20 }}>
                      <defs>
                        <linearGradient id="colorKn" x1="0" y1="0" x2="0" y2="1">
                          <stop offset="5%" stopColor="#4f46e5" stopOpacity={0.4}/>
                          <stop offset="95%" stopColor="#4f46e5" stopOpacity={0}/>
                        </linearGradient>
                      </defs>
                      <CartesianGrid strokeDasharray="3 3" stroke="#cbd5e1" opacity={0.6} />
                      <XAxis 
                        dataKey="webRegression" 
                        stroke="#64748b" 
                        label={{ value: 'Web Burned (mm)', position: 'insideBottom', offset: -10, fill: '#64748b' }} 
                        tick={{fontSize: 12}}
                      />
                      <YAxis 
                        yAxisId="left" 
                        stroke="#6366f1" 
                        label={{ value: 'Kn Ratio', angle: -90, position: 'insideLeft', fill: '#6366f1', offset: 0 }}
                        tick={{fontSize: 12}}
                      />
                      <YAxis 
                        yAxisId="right" 
                        orientation="right" 
                        stroke="#10b981" 
                        label={{ value: 'Web Remaining (mm)', angle: 90, position: 'insideRight', fill: '#10b981', offset: 0 }}
                        tick={{fontSize: 12}}
                      />
                      <Tooltip 
                        contentStyle={{ backgroundColor: '#ffffff', borderColor: '#e2e8f0', borderRadius: '8px', color: '#0f172a' }}
                        labelFormatter={(label) => `Web Burned: ${label} mm`}
                      />
                      <Legend verticalAlign="top" height={36}/>
                      <Area 
                        yAxisId="left"
                        type="monotone" 
                        dataKey="kn" 
                        name="Kn Ratio"
                        stroke="#4f46e5" 
                        strokeWidth={3} 
                        fillOpacity={1} 
                        fill="url(#colorKn)" 
                        animationDuration={1500}
                      />
                      <Line 
                        yAxisId="right"
                        type="monotone" 
                        dataKey="webRemaining" 
                        name="Web Remaining"
                        stroke="#10b981" 
                        strokeWidth={3} 
                        dot={false}
                        animationDuration={1500}
                      />
                    </ComposedChart>
                  </ResponsiveContainer>
                 ) : (
                   <div className="h-full flex items-center justify-center text-slate-500">Run evaluation to view data</div>
                 )}
              </div>
            </div>

            {/* Chart 4: Burn Rate Law */}
            <div className="bg-slate-100 border border-slate-200 rounded-xl p-6 shadow-sm h-[450px] flex flex-col">
              <div className="flex justify-between items-center mb-6">
                <div className="flex items-center space-x-2 text-orange-500">
                  <Flame className="w-5 h-5" />
                  <h3 className="text-lg font-bold text-slate-900">Burn Rate Law</h3>
                </div>
                {results && (
                   <div className="flex space-x-4 text-xs font-mono text-slate-500">
                      <span>r = {inputs.burnRateA} * P^{inputs.burnRateN}</span>
                   </div>
                )}
              </div>

              <div className="flex-1 w-full min-h-0">
                 {results ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <LineChart data={results.burnRateData} margin={{ top: 10, right: 30, left: 10, bottom: 20 }}>
                      <CartesianGrid strokeDasharray="3 3" stroke="#cbd5e1" opacity={0.6} />
                      <XAxis 
                        dataKey="pressure" 
                        type="number"
                        domain={[0, 10]}
                        stroke="#64748b" 
                        label={{ value: 'Pressure (MPa)', position: 'insideBottom', offset: -10, fill: '#64748b' }} 
                        tick={{fontSize: 12}}
                        tickCount={6}
                      />
                      <YAxis 
                        stroke="#f97316" 
                        label={{ value: 'Burn Rate (mm/s)', angle: -90, position: 'insideLeft', fill: '#f97316', offset: 0 }}
                        tick={{fontSize: 12}}
                      />
                      <Tooltip 
                        contentStyle={{ backgroundColor: '#ffffff', borderColor: '#e2e8f0', borderRadius: '8px', color: '#0f172a' }}
                        formatter={(value: number) => [value, 'mm/s']}
                        labelFormatter={(label) => `Pressure: ${label} MPa`}
                      />
                      <Line 
                        type="monotone" 
                        dataKey="rate" 
                        name="Burn Rate"
                        stroke="#f97316" 
                        strokeWidth={3} 
                        dot={false}
                        animationDuration={1500}
                      />
                    </LineChart>
                  </ResponsiveContainer>
                 ) : (
                   <div className="h-full flex items-center justify-center text-slate-500">Run evaluation to view data</div>
                 )}
              </div>
            </div>
          </div>
        )}
      </main>
    </div>
  );
};

export default App;
