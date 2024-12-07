extern crate nalgebra as na;
use minifb::{Key, Window, WindowOptions};
use na::{Rotation3, Unit, Vector3};
use std::f32::consts::PI;

struct Sphere {
    center: Vector3<f32>,
    radius: f32,
}

struct PhysicalProperties {
    ambient: Vector3<f32>,
    diffuse: Vector3<f32>,
    specular: Vector3<f32>,
}

struct PhysicalSphere {
    sphere: Sphere,
    physical_properties: PhysicalProperties,
    shininess: u32,
    reflection: f32,
}

struct PhysicalLight {
    position: Vector3<f32>,
    physical_properties: PhysicalProperties,
}

fn normalize(vector: Vector3<f32>) -> Vector3<f32> {
    vector / vector.norm()
}

fn nearest_intersected_object(
    objects: &Vec<PhysicalSphere>,
    ray_origin: Vector3<f32>,
    ray_direction: Vector3<f32>,
) -> (Option<&PhysicalSphere>, f32) {
    let mut min_distance = f32::INFINITY;
    let mut nearest_object = None;

    for object in objects {
        if let Some(distance) = sphere_intersect(
            object.sphere.center,
            object.sphere.radius,
            ray_origin,
            ray_direction,
        ) {
            if distance < min_distance {
                min_distance = distance;
                nearest_object = Some(object);
            }
        }
    }

    (nearest_object, min_distance)
}

fn sphere_intersect(
    center: Vector3<f32>,
    radius: f32,
    ray_origin: Vector3<f32>,
    ray_direction: Vector3<f32>,
) -> Option<f32> {
    let oc = ray_origin - center;
    let b = 2.0 * ray_direction.dot(&oc);
    let c = oc.dot(&oc) - radius * radius;
    let delta = b * b - 4.0 * c;

    if delta > 0.0 {
        let t1 = (-b + delta.sqrt()) / 2.0;
        let t2 = (-b - delta.sqrt()) / 2.0;
        if t1 > 0.0 && t2 > 0.0 {
            return Some(t1.min(t2));
        }
    }
    None
}

fn reflected(vector: Vector3<f32>, axis: Vector3<f32>) -> Vector3<f32> {
    vector - 2.0 * vector.dot(&axis) * axis
}

fn linspace(start: f32, stop: f32, num: usize) -> Vec<f32> {
    let step = (stop - start) / (num - 1) as f32;
    (0..num).map(|i| start + step * i as f32).collect()
}

fn main() {
    let width = 1280;
    let height = 720;

    let max_depth = 2;

    let mut camera_pos = Vector3::new(0.0, 0.0, 1.0);
    let mut camera_direction = Vector3::new(0.0, 0.0, -1.0).normalize();
    let mut camera_up = Vector3::new(0.0, 1.0, 0.0);
    let mut camera_right = camera_direction.cross(&camera_up).normalize();
    camera_up = camera_right.cross(&camera_direction).normalize(); // Ensure orthogonality

    let camera_speed = 0.1;
    let rotation_speed = camera_speed;

    let ratio = width as f32 / height as f32;

    let screen = [-1., 1. / ratio, 1., -1. / ratio];

    let light = PhysicalLight {
        position: Vector3::new(5., 5., 5.),
        physical_properties: PhysicalProperties {
            ambient: Vector3::new(1., 1., 1.),
            diffuse: Vector3::new(1., 1., 1.),
            specular: Vector3::new(1., 1., 1.),
        },
    };

    let objects = vec![
        PhysicalSphere {
            sphere: Sphere {
                center: Vector3::new(-0.2, 0.0, -1.0),
                radius: 0.7,
            },
            physical_properties: PhysicalProperties {
                ambient: Vector3::new(0.1, 0.0, 0.0),
                diffuse: Vector3::new(0.7, 0.0, 0.0),
                specular: Vector3::new(1.0, 1.0, 1.0),
            },
            shininess: 100,
            reflection: 0.5,
        },
        PhysicalSphere {
            sphere: Sphere {
                center: Vector3::new(0.1, -0.3, 0.0),
                radius: 0.1,
            },
            physical_properties: PhysicalProperties {
                ambient: Vector3::new(0.1, 0.0, 0.1),
                diffuse: Vector3::new(0.7, 0.0, 0.7),
                specular: Vector3::new(1.0, 1.0, 1.0),
            },
            shininess: 100,
            reflection: 0.5,
        },
        PhysicalSphere {
            sphere: Sphere {
                center: Vector3::new(-0.3, 0.0, 0.0),
                radius: 0.15,
            },
            physical_properties: PhysicalProperties {
                ambient: Vector3::new(0.0, 0.1, 0.0),
                diffuse: Vector3::new(0.0, 0.6, 0.0),
                specular: Vector3::new(1.0, 1.0, 1.0),
            },
            shininess: 100,
            reflection: 0.5,
        },
        PhysicalSphere {
            sphere: Sphere {
                center: Vector3::new(0.0, -9000.0, 0.0),
                radius: 9000.0 - 0.7,
            },
            physical_properties: PhysicalProperties {
                ambient: Vector3::new(0.1, 0.1, 0.1),
                diffuse: Vector3::new(0.6, 0.6, 0.6),
                specular: Vector3::new(1.0, 1.0, 1.0),
            },
            shininess: 40,
            reflection: 0.5,
        },
    ];

    let mut buffer: Vec<u32> = vec![0; width * height]; // 32-bit color buffer

    let mut window = Window::new(
        "Real-Time Ray Tracing",
        width,
        height,
        WindowOptions::default(),
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });

    while window.is_open() && !window.is_key_down(Key::Escape) {
        // Handle camera movement
        if window.is_key_down(Key::W) {
            camera_pos += camera_direction * camera_speed; // Move forward
        }
        if window.is_key_down(Key::S) {
            camera_pos -= camera_direction * camera_speed; // Move backward
        }
        if window.is_key_down(Key::A) {
            camera_pos -= camera_right * camera_speed; // Move left
        }
        if window.is_key_down(Key::D) {
            camera_pos += camera_right * camera_speed; // Move right
        }
        if window.is_key_down(Key::Space) {
            camera_pos += camera_up * camera_speed; // Move up
        }
        if window.is_key_down(Key::LeftShift) {
            camera_pos -= camera_up * camera_speed; // Move down
        }

        // Handle camera rotation
        let rotation_matrix_yaw =
            Rotation3::from_axis_angle(&Unit::new_normalize(camera_up), -rotation_speed);
        let rotation_matrix_pitch =
            Rotation3::from_axis_angle(&Unit::new_normalize(camera_right), rotation_speed);

        if window.is_key_down(Key::Left) {
            camera_direction = rotation_matrix_yaw.inverse() * camera_direction;
            camera_right = camera_direction.cross(&camera_up).normalize();
        }
        if window.is_key_down(Key::Right) {
            camera_direction = rotation_matrix_yaw * camera_direction;
            camera_right = camera_direction.cross(&camera_up).normalize();
        }
        if window.is_key_down(Key::Up) {
            camera_direction = rotation_matrix_pitch * camera_direction;
            camera_up = camera_right.cross(&camera_direction).normalize();
        }
        if window.is_key_down(Key::Down) {
            camera_direction = rotation_matrix_pitch.inverse() * camera_direction;
            camera_up = camera_right.cross(&camera_direction).normalize();
        }

        for (i, y) in linspace(screen[1], screen[3], height).iter().enumerate() {
            for (j, x) in linspace(screen[0], screen[2], width).iter().enumerate() {
                let pixel = camera_pos + camera_direction + camera_right * *x + camera_up * *y;
                let mut direction = (pixel - camera_pos).normalize();

                let mut origin = camera_pos;
                let mut color = Vector3::new(0.0, 0.0, 0.0);
                let mut reflection = 1.0;

                for _ in 0..max_depth {
                    let (nearest_object, min_distance) =
                        nearest_intersected_object(&objects, origin, direction);
                    if let Some(nearest_object) = nearest_object {
                        let intersection = origin + direction * min_distance;
                        let normal_to_surface =
                            (intersection - nearest_object.sphere.center).normalize();
                        let shifted_point = intersection + normal_to_surface * 1e-5;
                        let intersection_to_light = (light.position - shifted_point).normalize();

                        let (_, min_distance_to_light) = nearest_intersected_object(
                            &objects,
                            shifted_point,
                            intersection_to_light,
                        );
                        let intersection_to_light_distance = (light.position - intersection).norm();
                        let is_shadowed = min_distance_to_light < intersection_to_light_distance;

                        if is_shadowed {
                            break;
                        }

                        let mut illumination = Vector3::new(0.0, 0.0, 0.0);

                        // Ambient
                        illumination += nearest_object
                            .physical_properties
                            .ambient
                            .component_mul(&light.physical_properties.ambient);

                        // Diffuse
                        let diffuse_intensity =
                            intersection_to_light.dot(&normal_to_surface).max(0.0);
                        illumination += nearest_object
                            .physical_properties
                            .diffuse
                            .component_mul(&light.physical_properties.diffuse)
                            * diffuse_intensity;

                        // Specular
                        let intersection_to_camera = (camera_pos - intersection).normalize();
                        let half_vector =
                            (intersection_to_light + intersection_to_camera).normalize();
                        let specular_intensity = normal_to_surface
                            .dot(&half_vector)
                            .max(0.0)
                            .powf(nearest_object.shininess as f32 / 4.0);
                        illumination += nearest_object
                            .physical_properties
                            .specular
                            .component_mul(&light.physical_properties.specular)
                            * specular_intensity;

                        // Reflection
                        color += illumination * reflection;
                        reflection *= nearest_object.reflection;

                        origin = shifted_point;
                        direction = reflected(direction, normal_to_surface);
                    } else {
                        break;
                    }
                }

                // Set the pixel color in the image, clamping the values
                let r = (color.x * 255.0).min(255.0).max(0.0) as u8;
                let g = (color.y * 255.0).min(255.0).max(0.0) as u8;
                let b = (color.z * 255.0).min(255.0).max(0.0) as u8;
                let color: u32 = (255 << 24) | ((r as u32) << 16) | ((g as u32) << 8) | (b as u32);

                // Write the color to the buffer
                buffer[i * width + j] = color;
            }
        }
        window.update_with_buffer(&buffer, width, height).unwrap();
    }
}
