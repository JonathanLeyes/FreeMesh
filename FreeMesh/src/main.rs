use std::io::LineWriter;
//#[macro_use]
//extern crate nalgebra as na;
use std::ops;
use std::fmt;
use plotters::data;
use plotters::prelude::*;

use std::fs::OpenOptions;
use std::io::Write;

use std::time::{Duration, Instant};
use std::sync::{Arc, Mutex};


struct Point {
    x: f64,
    y: f64,
}

struct Scalar { value: f64}

// Implement basic math for scalars
impl ops::Add<Scalar> for Scalar {
    type Output = Self;
    fn add(self, rhs: Scalar) -> Self::Output {
        Self {
            value: self.value+rhs.value
        }
    }
}

impl ops::Sub<Scalar> for Scalar {
    type Output = Self;
    fn sub(self, rhs: Scalar) -> Self::Output {
        Self {
            value: self.value-rhs.value
        }
    }
}

impl ops::Mul<Scalar> for Scalar {
    type Output = Self;
    fn mul(self, rhs: Scalar) -> Self::Output {
        Self {
            value: self.value*rhs.value
        }
    }
}

impl ops::Div<Scalar> for Scalar {
    type Output = Self;
    fn div(self, rhs: Scalar) -> Self::Output {
        if rhs.value == 0.0 {
            panic!("Cannot divide by zero!")
        }
        Self {
            value: self.value/rhs.value
        }
    }
}

impl ops::Mul<Point> for Scalar {
    type Output = Point;
    fn mul(self, rhs: Point) -> Self::Output {
        Point {
            x: self.value * rhs.x,
            y: self.value * rhs.y,
        }
    }
}

// Implement display of scalar
impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

// Implement copy of scalar
impl Copy for Scalar {}

impl Clone for Scalar {
    fn clone(&self) -> Scalar {
        *self
    }
}

impl Point {
    fn new(x: f64, y: f64) -> Point {
        Point {x:x, y:y}
    }
    fn subtract(P1: Point, P2: Point) -> Point {
        Point{x:P1.x-P2.x, y:P1.y-P2.y}
    }
    fn abs(P: Point) -> Scalar {
        Scalar { value: (P.x*P.x + P.y*P.y).sqrt()}
    }
    fn scale(scaler: f64, P:Point) -> Point {
        Point {x:P.x*scaler, y:P.y*scaler}
    }
}

// Implement basic math for points
impl ops::Add for Point {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl ops::Sub for Point {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl ops::Mul<Scalar> for Point {
    type Output = Self;
    fn mul(self, rhs: Scalar) -> Self::Output {
        Self {
            x: self.x * rhs.value,
            y: self.y * rhs.value,
        }
    }
}

impl ops::Mul<Point> for Point {
    type Output = Self;
    fn mul(self, rhs: Point) -> Self::Output {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
        }
    }
}

// Implement display of points
impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

// Implement copy of points
impl Copy for Point {}

impl Clone for Point {
    fn clone(&self) -> Point {
        *self
    }
}

fn bond_vector(P1:Point, P2:Point) -> f64 {
    return f64::sqrt(f64::powf(P1.x - P2.x, 2.) + f64::powf(P1.y - P2.y, 2.));
}

fn bond_vec_arr(P1:Point, P:&Vec<Point>, id:usize) -> Vec<f64>{
    //Bond vector:  ξ = X`-X

    // input: Location of Point and Locations of Pointcollection
    // return: Collection of distances

    let mut P_coll = vec![0.0; P.len()];
    let mut c = 0;
    for i in P {
        if id != c {
            P_coll[c] = f64::sqrt(f64::powf(i.x - P1.x, 2.) + f64::powf(i.y - P1.y, 2.));
        }
        else {
            P_coll[c] = 0.0
        }
        
        c += 1;
    }
    return P_coll;
}

fn bond_displacement(P1u:f64, P2u:f64) -> f64 {
    return P2u - P1u;
}

fn bond_displacement_arr(P1u:f64, Pu:&Vec<f64>, id:usize) -> Vec<f64> {
    // bond displacement: η(t) = u(X`,t) - u(X,t)

    // input: Displacment of Point and Displacement of Pointcollection
    // return: Collection of displacement-differences

    let mut P_coll = vec![0.0; Pu.len()];
    let mut c = 0;
    for i in Pu {
        if id != c {
            P_coll[c] = i-P1u;
        }
        else {
            P_coll[c] = 0.;
        }
        
        c += 1;
    }

    return P_coll;
}

fn main() {
    // Coordinates of Points
    let start = Instant::now();
    let domain = [0., 10., 0., 10.];
    const nx:usize = 80;
    const ny:usize = 80;
    const node_count:usize = nx*ny;

    let dx = domain[1]/nx as f64;
    let dy = domain[3]/ny as f64;

    let delta = 3.0;
    let E = 210.;
    let nu = 0.3;
    let K = E/(3.-6.*nu);


    let mut p_col = vec![Point{x:0., y:0.}; node_count];
    let mut p_vol_col = vec![1.0; node_count];
    let mut p_u_col = vec![0.; node_count];

    p_u_col[300] = 0.5; // definie an initial local displacement

    // populate domain with nodes
    let mut c = 0;

    for i in 0..nx as i32 {
        for j in 0..ny as i32 {
            p_col[c] = Point {x:domain[0]+dx*i as f64, y:dy*j as f64};
            c += 1;
        }
    }


    // Energy density function W(X,t) + 1/2 ∫ ω(X`-X)w[u(X`, t) - u(X,t), X`-X] dV_X`
    // micropotential: w[u(X`, t) - u(X,t)]
    // Weight function: ω(X`-X)
    // Bond vector:  ξ = X`-X
    // bond displacement: η(t) = u(X`,t) - u(X,t)
    // bond distance vector: r(t) = η(t) + ξ

    // Discrete form: W_i = Σ ω(ξ_ij) V_j w_ij(r_ij, ω_ij)



    // Linear Elasticity

    // f(s,ξ) = -cs/ξ, c: proportionality constant
    // micropotential: w(s) = -∫f(s,ξ)dη = 1/2cs^2
    // weight function: ω(ξ_ij) = 1 if ξ_ij <= δ || ω(ξ_ij) = 0 if ξ_ij > δ
    // c=6K/(π δ^3), K: bulk modulus K=E/(3-6 \nu)
    // oder c_i = 18K/(Σ V_j), j ∈ H_δ

    let mut W = Arc::new(Mutex::new(vec![0.;node_count]));
    let mut W_max = 0.0;
    let mut W_min = 1e+16;

    let mut file = Arc::new(Mutex::new(OpenOptions::new().create_new(true).append(true).open("raw_data.txt").expect("cannot open.")));

    
    //let mut datarows = vec![];

    // main loop
    use rayon::prelude::*;

    p_col.par_iter().enumerate().for_each(|(counter, i)| {
        let xi = bond_vec_arr(*i, &p_col, counter);
        let eta = bond_displacement_arr(p_u_col[counter], &p_u_col, counter);

        let omega: Vec<bool> = xi.iter().map(|x|x <= &delta).collect();
        let mut v_sum = 0.;
        for (c, v) in p_vol_col.iter().enumerate() {
            if counter != c {
                if omega[c] == true {
                    v_sum += *v;
                }
            }
        }
        let c = 6.*K/(3.14*f64::powf(delta, 3.));
        for j in 0..node_count
        {
            if counter != j {
                let s = f64::abs(xi[j] + eta[j])/f64::abs(xi[j]);
                if xi[j] <= delta
                    {
                        W.lock().unwrap()[counter] += p_vol_col[counter]*0.5*c*f64::powf(s, 2.);    // energy density trasnfered to linear elasticity
                        /*if W.lock().unwrap()[counter] > W_max {
                            //W_max = W.lock().unwrap()[counter];
                        }
                        if W.lock().unwrap()[counter] < W_min {
                            //W_min = W.lock().unwrap()[counter];
                        }*/
                        
                    }
            }
        }
        //datarows.push((i.x, i.y, W.lock().unwrap()[counter]));

        let mut str = format!("{}, {}, {}\n", i.x, i.y, W.lock().unwrap()[counter]);
            file.lock().unwrap().write_all( str.as_bytes()).expect("write failed");
    });
/* 
    for (counter, i) in p_col.iter().enumerate() {
        
        let xi = bond_vec_arr(*i, &p_col, counter);
        let eta = bond_displacement_arr(p_u_col[counter], &p_u_col, counter);

        let omega: Vec<bool> = xi.iter().map(|x|x <= &delta).collect();
        let mut v_sum = 0.;
        for (c, v) in p_vol_col.iter().enumerate() {
            if counter != c {
                if omega[c] == true {
                    v_sum += *v;
                }
            }
        }

        //let c = 18.*K/v_sum;
        let c = 6.*K/(3.14*f64::powf(delta, 3.));
        

        for j in 0..node_count
        {
            if counter != j {
                let s = f64::abs(xi[j] + eta[j])/f64::abs(xi[j]);
                if xi[j] <= delta
                    {
                        W[counter] += p_vol_col[counter]*0.5*c*f64::powf(s, 2.);    // energy density trasnfered to linear elasticity
                        if W[counter] > W_max {
                            W_max = W[counter];
                        }
                        if W[counter] < W_min {
                            W_min = W[counter];
                        }
                        
                    }
            }
        }
        datarows.push((i.x, i.y, W[counter]));

        let mut str = format!("{}, {}, {}\n", i.x, i.y, W[counter]);
            file.write_all( str.as_bytes()).expect("write failed");
        
    }
*/
    let duration = start.elapsed();
    println!("Time elapsed: {:?}", duration);


/*
    println!("W1: {}", W[0]);
    println!("W2: {}", W[1]);
    println!("W3: {}", W[2]);
    println!("W4: {}", W[3]);
    println!("W5: {}", W[4]);
    println!("W100: {}", W[99]);

*/

    let P1_coor = Point{x:0., y:0.};
    let P1_mass = Scalar {value: 10.};
    let P1_rho:Scalar = Scalar {value: 4.};
    let P1_u:Scalar = Scalar {value: 10.};

    let P2_coor = Point{x:0., y:1.};
    let P2_mass = Scalar {value: 5.};
    let P2_rho = Scalar {value: 4.};
    let P2_u:Scalar = Scalar {value: 8.};

    let P1_u_new = Point::abs(P1_coor-P1_coor)*(P1_mass/P1_rho*P1_u) + Point::abs(P1_coor-P2_coor)*(P2_mass/P2_rho*P2_u); // definitionof u(x,t_n) = SUM_i ( m_i * u_i^n/rho_i (|x-x_i|))
    let P2_u_new = Point::abs(P2_coor-P1_coor)*(P1_mass/P1_rho*P1_u) + Point::abs(P2_coor-P2_coor)*(P2_mass/P2_rho*P2_u);
    //let P2_u_new = P1_mass*P1_u/P1_rho*(na::Vector2::abs(P1_cor-P1_cor)) + P2_mass*P2_u/P2_rho*(na::Vector2::abs(P2_cor-P2_cor));


    println!("P1_u: {} -> {}", P1_u, P1_u_new);
    println!("P2_u: {} -> {}", P2_u, P2_u_new);
}