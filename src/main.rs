#![feature(core)]
extern crate nalgebra as na;
extern crate term;
extern crate csv;
extern crate rustc_serialize;
use std::io::prelude::*;
use na::{ Mat2, Vec2 };
use std::cmp::{ partial_min, partial_max };
use std::f64::consts;

fn build_region(lat: f64, lon: f64, h: f64, bear: f64, x: i64, y: i64) -> (Box<Fn(f64,f64) -> bool>, f64, f64, f64, f64) {
    const RESOLUTION_X : f64 = 20./360.0*consts::PI_2;
    const RESOLUTION_Y : f64 = 5./360.0*consts::PI_2;
    let long_units = h*RESOLUTION_X.tan()/8.0/110574.0;
    let lat_units = h*RESOLUTION_Y.tan()/2.0/110574.0;
    //println!("{} {}", long_units, lat_units);
    let rotation = Mat2::new(bear.cos(), -bear.sin(), bear.sin(), bear.cos());
    let tl = rotation*Vec2::new(((x-8) as f64)*long_units, ((2-y) as f64)*lat_units) + Vec2::new(lon, lat);
    let tr = rotation*Vec2::new(((x-8+1) as f64)*long_units, ((2-y) as f64)*lat_units) + Vec2::new(lon, lat); 
    let bl = rotation*Vec2::new(((x-8) as f64)*long_units, ((2-y-1) as f64)*lat_units) + Vec2::new(lon, lat); 
    let br = rotation*Vec2::new(((x-8+1) as f64)*long_units, ((2-y-1) as f64)*lat_units) + Vec2::new(lon, lat);

    //println!("{}, {} {}, {} {}, {} {}, {}", tl.x, tl.y, tr.x, tr.y, br.x, br.y, bl.x, bl.y);
    
    let max_x = partial_max(partial_max(tl.x, br.x).unwrap(),partial_max(bl.x,tr.x).unwrap()).unwrap();
    let max_y = partial_max(partial_max(tl.y, br.y).unwrap(),partial_max(bl.y,tr.y).unwrap()).unwrap();
    let min_x = partial_min(partial_min(tl.x, br.x).unwrap(),partial_min(bl.x,tr.x).unwrap()).unwrap();
    let min_y = partial_min(partial_min(tl.y, br.y).unwrap(),partial_min(bl.y,tr.y).unwrap()).unwrap();

    /*fn rect_area(x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) -> f64 {
        let a = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
        let b = ((x3 - x2).powi(2) + (y3 - y2).powi(2)).sqrt();
        let area = a*b;
        area
    }*/
    let rect_are = (long_units*lat_units).abs();
    fn tring_area(x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) -> f64 {
        let a = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
        let b = ((x3 - x2).powi(2) + (y3 - y2).powi(2)).sqrt();
        let c = ((x3 - x1).powi(2) + (y3 - y1).powi(2)).sqrt();

        let s = (a + b + c)/2.;
        let area = (s*(s-a)*(s-b)*(s-c)).sqrt();
        area
    }
    

    (Box::new(move |x,y| { 
        //println!("{} {}", rect_area(tl.x,tl.y,tr.x,tr.y,br.x,br.y), rect_are);
        if (tring_area(x, y, tl.x, tl.y, bl.x, bl.y) + tring_area(x, y, tl.x, tl.y, tr.x, tr.y) + tring_area(x, y, tr.x, tr.y, br.x, br.y) + tring_area(x, y, br.x, br.y, bl.x, bl.y) - rect_are).abs() < 0.000000000000000001 {
            true 
        } else { 
            false 
        }
    }), max_x, min_x, max_y, min_y)
}

fn build_all_regions(temps: Vec<Vec<Vec<i64>>>, env_data: Vec<(f64, f64, f64, f64)>) -> (Vec<Vec<(Box<Fn(f64,f64) -> bool>,i64)>>, f64,f64,f64,f64) {
   let mut max_lat = std::f64::MIN;
   let mut min_lat = std::f64::MAX;
   let mut max_long = std::f64::MIN;
   let mut min_long = std::f64::MAX;
   let mut h_count = 0usize;
   let mut shot: Vec<(Box<Fn(f64,f64) -> bool>,i64)> = vec![];
   let mut shots: Vec<Vec<(Box<Fn(f64,f64) -> bool>,i64)>> = vec![];
    for h in &env_data {
        for i in 0..4 {
            for j in 0..16 {
                let (a,b,c,d,e) = build_region(h.0,h.1, h.2, h.3, j, i);
                shot.push((a, temps[h_count][i as usize][j as usize]));
                max_lat = partial_max(max_lat, d).unwrap();
                min_lat = partial_min(min_lat, e).unwrap();
                max_long = partial_max(max_long, b).unwrap();
                min_long = partial_min(min_long, c).unwrap();
            }
        }
        //println!("{},{},{},{}", max_lat, min_lat, max_long, min_long);
        shots.push(shot);
        shot = vec![];
        h_count += 1;
    }
    (shots, max_long, min_long, max_lat, min_lat)
}

struct Region {
    temp: i64,
    mid_x: f64,
    mid_y: f64,
}

#[derive(RustcDecodable)]
struct Record {
    t11: i64, t12: i64, t13: i64, t14: i64, t15: i64, t16: i64, t17: i64, t18: i64, t19: i64, t110: i64, t111: i64, t112: i64, t113: i64, t114: i64, t115: i64, t116: i64,
    t21: i64, t22: i64, t23: i64, t24: i64, t25: i64, t26: i64, t27: i64, t28: i64, t29: i64, t210: i64, t211: i64, t212: i64, t213: i64, t214: i64, t215: i64, t216: i64,
    t31: i64, t32: i64, t33: i64, t34: i64, t35: i64, t36: i64, t37: i64, t38: i64, t39: i64, t310: i64, t311: i64, t312: i64, t313: i64, t314: i64, t315: i64, t316: i64,
    t41: i64, t42: i64, t43: i64, t44: i64, t45: i64, t46: i64, t47: i64, t48: i64, t49: i64, t410: i64, t411: i64, t412: i64, t413: i64, t414: i64, t415: i64, t416: i64,
    bear: f64, h: f64, lat: f64, long: f64,
}

fn get_data(path: String) -> (Vec<Vec<Vec<i64>>>, Vec<(f64,f64,f64,f64)>) {
     let mut temps: Vec<Vec<Vec<i64>>> = vec![];
   let mut env_data: Vec<(f64,f64,f64,f64)> = vec![];
   let mut rdr = csv::Reader::from_file(path).unwrap();

    for record in rdr.decode() {
         let record: Record = record.unwrap();
       temps.push(vec![vec![record.t11,record.t12,record.t13,record.t14,record.t15,record.t16,record.t17,record.t18,record.t19,record.t110,record.t111,record.t112,record.t113,record.t114,record.t115,record.t116],
                        vec![record.t21,record.t22,record.t23,record.t24,record.t25,record.t26,record.t27,record.t28,record.t29,record.t210,record.t211,record.t212,record.t213,record.t214,record.t215,record.t216],
                        vec![record.t31,record.t32,record.t33,record.t34,record.t35,record.t36,record.t37,record.t38,record.t39,record.t310,record.t311,record.t312,record.t313,record.t314,record.t315,record.t316],
                        vec![record.t41,record.t42,record.t43,record.t44,record.t45,record.t46,record.t47,record.t48,record.t49,record.t410,record.t411,record.t412,record.t413,record.t414,record.t415,record.t416]]);

        env_data.push((record.lat, record.long, record.h, (record.bear + 90.)/360.0*consts::PI_2));
    }
    (temps, env_data)
}

fn main() {
   const RESOLUTION: f64 = 100.0;
   const RESOLUTION_2: f64 = RESOLUTION*2.0;

   let (temps, env_data) = get_data("./my_data.csv".to_string());

   let (shots, max_long, min_long, max_lat, min_lat) = build_all_regions(temps, env_data);
    let mut map: Vec<Region> = vec![];
    let min_across = partial_min((max_long-min_long), (max_lat-min_lat)).unwrap();
    for i in 0..(((max_lat - min_lat)*RESOLUTION/min_across).round() as i64) {
         for j in 0..(((max_long - min_long)*RESOLUTION/min_across).round() as i64) {
             map.push(Region{ temp: std::i64::MAX, mid_x: min_long + min_across/RESOLUTION_2 + (j as f64)*min_across/RESOLUTION, mid_y: min_lat + min_across/RESOLUTION_2 + (i as f64)*min_across/RESOLUTION });
         }
     }

    for h in map.iter_mut() {
        for i in &shots {
            for j in i {
                if j.0(h.mid_x, h.mid_y) {
                    if h.temp == std::i64::MAX { h.temp = j.1; } else { h.temp = (j.1 + h.temp)/2;}
                }
            }
        }
    }

    let mut t = term::stdout().unwrap();
    for i in 0..(((max_lat - min_lat)*RESOLUTION/min_across).round() as usize) {
        for j in (i*(((max_long - min_long)*RESOLUTION/min_across).round() as usize))..((i+1usize)*(((max_long - min_long)*RESOLUTION/min_across).round() as usize)) {
            match map[j].temp {
                1...20 => t.fg(term::color::CYAN).unwrap(),
                -3...0 => t.fg(term::color::RED).unwrap(),
                -7...-4 => t.fg(term::color::MAGENTA).unwrap(),
                -11...-8 => t.fg(term::color::YELLOW).unwrap(),
                -15...-12 => t.fg(term::color::WHITE).unwrap(),
                _ => t.fg(term::color::BLACK).unwrap(),
            };
             (write!(t, ".")).unwrap();

        }
        println!("");
    }
}






