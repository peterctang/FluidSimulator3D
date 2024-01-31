/*#include "fluid.h"

Fluid::Fluid()
{
	diffusion = 0.00001f;
	viscosity = 0.000000f;
	buoyancy = 4.0f;
	vc_eps = 5.0f;
	_dt = DT;
	_isLightSelected = false;
	_isRendering = true;
	_isDrawSliceOutline = false;

	for (int i = 0; i < 10; i++)
		ClearBuffer(_buffers[i]);

	int i = 0;
	_density = _buffers[i++]; _densityTmp = _buffers[i++];
	_velX = _buffers[i++]; _velXTmp = _buffers[i++];
	_velY = _buffers[i++]; _velYTmp = _buffers[i++];
	_velZ = _buffers[i++]; _velZTmp = _buffers[i++];

	ClearSources();

	_lightPos[0] = -1.2f;
	_lightPos[1] = 0.2f;
	_lightPos[2] = 1.2f;
	_renderer = new Renderer(_density, RES);
	_renderer->SetLightPostion(_lightPos);
}

const float* Fluid::GetDensity()
{
	return _density;
}

bool Fluid::LightSelected(double mouseX, double mouseY)
{
	GLdouble mvMatrix[16], projMatrix[16];
	GLint viewportMatrix[4];
	double objX = 0, objY = 0, objZ = 0;
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewportMatrix);

	//point on 'near' plane
	gluUnProject(mouseX, _winY - mouseY, 0.0, mvMatrix, projMatrix,
		viewportMatrix, &objX, &objY, &objZ);
	Eigen::Vector3f ptNear(objX, objY, objZ);

	//point on 'far' plane
	gluUnProject(mouseX, _winY - mouseY, 1.0, mvMatrix, projMatrix,
		viewportMatrix, &objX, &objY, &objZ);
	Eigen::Vector3f ptFar(objX, objY, objZ);

#if 1
	//calculate distance between point and line
	Eigen::Vector3f crossProduct = (_lightPos - ptNear).cross(_lightPos - ptFar);
	float dist = crossProduct.norm() /
		(ptFar - ptNear).norm();
#else
	//another method: calculate distance between point and line
	Eigen::Vector3f lineDir = (ptFar - ptNear);
	lineDir.normalize();
	Eigen::Vector3f pointDir = (_lightPos - ptNear);
	float t = pointDir.dot(lineDir);
	Eigen::Vector3f projection = ptNear + (lineDir * t);
	float dist = (projection - _lightPos).norm();
#endif
	if (dist < 0.1)
		return true;

	return false;
}

void Fluid::MouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double mouseX, mouseY;
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		::glfwGetCursorPos(window, &mouseX, &mouseY);
		if (action == GLFW_PRESS) {
			_isLeftKeyPressed = true;

			if (LightSelected(mouseX, mouseY)) {
				std::cout << "light selected" << std::endl;
				_isLightSelected = true;
				_arcball.StartDragging(mouseX, mouseY);
			}
			_arcball.StartRotation(mouseX, mouseY);
		}
		else if (action == GLFW_RELEASE) {
			_isLeftKeyPressed = false;
			_isLightSelected = false;
			_arcball.StopRotation();
			_arcball.StopDragging();
		}
	}
	else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
		_isMiddleKeyPressed = (action == GLFW_PRESS);
	}
	else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
#if 0
		if (action == GLFW_PRESS) {
			_isRightKeyPressed = true;
			::glfwGetCursorPos(window, &mouseX, &mouseY);
			_arcball.StartZooming(mouseX, mouseY);
		}
		else if (action == GLFW_RELEASE) {
			_isRightKeyPressed = false;
			_arcball.StopZooming();
		}
#endif
	}
}

void Fluid::MouseMotion(GLFWwindow* window, double nx, double ny)
{
	if (_isLeftKeyPressed && _isCtrlPressed) {
	}
	else if (_isLeftKeyPressed && _isLightSelected) {
		_lightPos += _arcball.UpdateDragging(nx, ny);
		_renderer->SetLightPostion(_lightPos);
	}
	else if (_isLeftKeyPressed) {
		_arcball.UpdateRotation(nx, ny);
	}
#if 0
	else if (_isRightKeyPressed) {
		_arcball.UpdateZooming(nx, ny);
	}
#endif
}

void Fluid::Keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch (key) {
		case GLFW_KEY_S:		// Rendering on/off
			_isRendering = !_isRendering;
			_renderer->SetRendering(_isRendering);
			break;
		case GLFW_KEY_W:		//wireframe on/off
			_isDrawSliceOutline = !_isDrawSliceOutline;
			_renderer->SetSliceOutline(_isDrawSliceOutline);
			break;
		}
	}
}

Fluid::~Fluid()
{
}

#define N 38
#define BOUNDARY
#undef BOUNDARY
void Fluid::set_bnd(int N, int b, float* x) {
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			x[IX(0, i, j)] = (b == 1) ? -x[IX(1, i, j)] : x[IX(1, i, j)];
			x[IX(N + 1, i, j)] = (b == 1) ? -x[IX(N, i, j)] : x[IX(N, i, j)];
			x[IX(i, 0, j)] = (b == 2) ? -x[IX(i, 1, j)] : x[IX(i, 1, j)];
			x[IX(i, N + 1, j)] = (b == 2) ? -x[IX(i, N, j)] : x[IX(i, N, j)];
			x[IX(i, j, 0)] = (b == 3) ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
			x[IX(i, j, N + 1)] = (b == 3) ? -x[IX(i, j, N)] : x[IX(i, j, N)];
		}
	}

	// Corners
	x[IX(0, 0, 0)] = 0.5 * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
	x[IX(0, 0, N + 1)] = 0.5 * (x[IX(1, 0, N + 1)] + x[IX(0, 1, N + 1)] + x[IX(0, 0, N)]);
	x[IX(0, N + 1, 0)] = 0.5 * (x[IX(1, N + 1, 0)] + x[IX(0, N, 0)] + x[IX(0, N + 1, 1)]);
	x[IX(0, N + 1, N + 1)] = 0.5 * (x[IX(1, N + 1, N + 1)] + x[IX(0, N, N + 1)] + x[IX(0, N + 1, N)]);
	x[IX(N + 1, 0, 0)] = 0.5 * (x[IX(N, 0, 0)] + x[IX(N + 1, 1, 0)] + x[IX(N + 1, 0, 1)]);
	x[IX(N + 1, 0, N + 1)] = 0.5 * (x[IX(N, 0, N + 1)] + x[IX(N + 1, 1, N + 1)] + x[IX(N + 1, 0, N)]);
	x[IX(N + 1, N + 1, 0)] = 0.5 * (x[IX(N, N + 1, 0)] + x[IX(N + 1, N, 0)] + x[IX(N + 1, N + 1, 1)]);
	x[IX(N + 1, N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1, N + 1)] + x[IX(N + 1, N, N + 1)] + x[IX(N + 1, N + 1, N)]);
}

void Fluid::add_source(int N, float* x, float* s, float dt)
{
	int i, size = (N + 2) * (N + 2) * (N + 2);
	for (i = 0; i < size; i++) x[i] += dt * s[i];
}

void Fluid::AddBuoyancy()
{
	int i;

	for (i = 0; i < SIZE; i++)
		_velY[i] += _density[i] * buoyancy * _dt;//FIXME
}

inline void Fluid::diffuse(int N, int b, float* x, float* x0, float diff, float dt)
{
	int i, j, k;
	float a = dt * diff * N * N;

	for (int n = 0; n < 20; n++) {
		for (i = 1; i <= N; i++) {
			for (j = 1; j <= N; j++) {
				for (k = 1; k <= N; k++) {
					x[IX(i, j, k)] = (x0[IX(i, j, k)] + a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] +
						x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] + x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / (1 + 6 * a);
				}
			}
		}
	}

	set_bnd(N, b, x);
}

inline void Fluid::advect(int N, int b, float* d, float* d0, float* u, float* v, float* w, float dt)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, r0, s1, t1, r1, dt0;

	dt0 = dt * N;

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			for (k = 1; k <= N; k++) {
				x = i - dt0 * u[IX(i, j, k)];
				y = j - dt0 * v[IX(i, j, k)];
				z = k - dt0 * w[IX(i, j, k)];

				// Clamp the coordinates
				if (x < 0.5) x = 0.5;
				if (x > N + 0.5) x = N + 0.5;
				i0 = (int)x;
				i1 = i0 + 1;

				if (y < 0.5) y = 0.5;
				if (y > N + 0.5) y = N + 0.5;
				j0 = (int)y;
				j1 = j0 + 1;

				if (z < 0.5) z = 0.5;
				if (z > N + 0.5) z = N + 0.5;
				k0 = (int)z;
				k1 = k0 + 1;

				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;
				r1 = z - k0;
				r0 = 1 - r1;

				d[IX(i, j, k)] = s0 * (t0 * (r0 * d0[IX(i0, j0, k0)] + r1 * d0[IX(i0, j0, k1)]) +
					t1 * (r0 * d0[IX(i0, j1, k0)] + r1 * d0[IX(i0, j1, k1)])) +
					s1 * (t0 * (r0 * d0[IX(i1, j0, k0)] + r1 * d0[IX(i1, j0, k1)]) +
						t1 * (r0 * d0[IX(i1, j1, k0)] + r1 * d0[IX(i1, j1, k1)]));
			}
		}
	}

	set_bnd(N, b, d);
}

void Fluid::project()
{
	int i, j, k;
	float h;
	h = 1.0 / N;
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			for (k = 1; k <= N; k++) {
				div[IX(i, j, k)] = -0.5 * h * (u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)] +
					v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)] +
					w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]);
				p[IX(i, j, k)] = 0;
			}
		}
	}
	set_bnd(N, 0, div);
	set_bnd(N, 0, p);
	for (int iter = 0; iter < 20; iter++) {
		for (i = 1; i <= N; i++) {
			for (j = 1; j <= N; j++) {
				for (k = 1; k <= N; k++) {
					p[IX(i, j, k)] = (div[IX(i, j, k)] + p[IX(i - 1, j, k)] + p[IX(i + 1, j, k)] +
						p[IX(i, j - 1, k)] + p[IX(i, j + 1, k)] + p[IX(i, j, k - 1)] + p[IX(i, j, k + 1)]) / 6;
				}
			}
		}
		set_bnd(N, 0, p);
	}
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			for (k = 1; k <= N; k++) {
				u[IX(i, j, k)] -= 0.5 * (p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]) / h;
				v[IX(i, j, k)] -= 0.5 * (p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]) / h;
				w[IX(i, j, k)] -= 0.5 * (p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]) / h;
			}
		}
	}
	set_bnd(N, 1, u);
	set_bnd(N, 2, v);
	set_bnd(N, 3, w);
}

void Fluid::VorticityConfinement()
{
	int ijk;
	//temp buffers
	float* curlX = _velXTmp, * curlY = _velYTmp, * curlZ = _velZTmp, * curl = _densityTmp;
	float dt0 = _dt * vc_eps;

	FOR_ALL_CELL{
		ijk = _I(i,j,k);
	// curlx = dw/dy - dv/dz
	curlX[ijk] = (_velZ[_I(i,j + 1,k)] - _velZ[_I(i,j - 1,k)]) * 0.5f -
		(_velY[_I(i,j,k + 1)] - _velY[_I(i,j,k - 1)]) * 0.5f;

	// curly = du/dz - dw/dx
	curlY[ijk] = (_velX[_I(i,j,k + 1)] - _velX[_I(i,j,k - 1)]) * 0.5f -
		(_velZ[_I(i + 1,j,k)] - _velZ[_I(i - 1,j,k)]) * 0.5f;

	// curlz = dv/dx - du/dy
	curlZ[ijk] = (_velY[_I(i + 1,j,k)] - _velY[_I(i - 1,j,k)]) * 0.5f -
		(_velX[_I(i,j + 1,k)] - _velX[_I(i,j - 1,k)]) * 0.5f;

	// curl = |curl|
	curl[ijk] = sqrtf(curlX[ijk] * curlX[ijk] +
			curlY[ijk] * curlY[ijk] +
			curlZ[ijk] * curlZ[ijk]);
	} END_FOR

		FOR_ALL_CELL{
			ijk = _I(i,j,k);
			float nX = (curl[_I(i + 1,j,k)] - curl[_I(i - 1,j,k)]) * 0.5f;
			float nY = (curl[_I(i,j + 1,k)] - curl[_I(i,j - 1,k)]) * 0.5f;
			float nZ = (curl[_I(i,j,k + 1)] - curl[_I(i,j,k - 1)]) * 0.5f;
			float len1 = 1.0f / (sqrtf(nX * nX + nY * nY + nZ * nZ) + 0.0000001f);
			nX *= len1;
			nY *= len1;
			nZ *= len1;
			_velX[ijk] += (nY * curlZ[ijk] - nZ * curlY[ijk]) * dt0;
			_velY[ijk] += (nZ * curlX[ijk] - nX * curlZ[ijk]) * dt0;
			_velZ[ijk] += (nX * curlY[ijk] - nY * curlX[ijk]) * dt0;
	} END_FOR
}

#define DIFFUSE
#define ADVECT

void Fluid::vel_step(int N, float* u, float* v, float* w, float* u0, float* v0, float* w0, float visc, float dt)
{
	// Add external forces (if any)
	add_source(N, u, u0, dt);
	add_source(N, v, v0, dt);
	add_source(N, w, w0, dt);

	// Swap velocity fields
	SWAP(u0, u);
	SWAP(v0, v);
	SWAP(w0, w);

	// Diffuse velocities
	diffuse(N, 1, u, u0, visc, dt);
	diffuse(N, 2, v, v0, visc, dt);
	diffuse(N, 3, w, w0, visc, dt);

	// Project the velocities to make it divergence-free
	project(N, u, v, w, u0, v0);

	// Swap velocity fields back
	SWAP(u0, u);
	SWAP(v0, v);
	SWAP(w0, w);

	// Advect velocities
	advect(N, 1, u, u0, u0, v0, w0, dt);
	advect(N, 2, v, v0, u0, v0, w0, dt);
	advect(N, 3, w, w0, u0, v0, w0, dt);

	// Project the velocities again
	project(N, u, v, w, u0, v0);
}

void Fluid::dens_step(int N, float* x, float* x0, float* u, float* v, float* w, float diff, float dt)
{
	add_source(N, x, x0, dt);
	SWAP(x0, x);
	diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x);
	advect(N, 0, x, x0, u, v, w, dt);
}

void Fluid::GenerateSmoke()
{
	const int centerY = RES / 4;
	const int centerZ = RES / 2;
	float dens = (rand() % 1000) / 1000.0f;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			if (hypot(i - centerY, j - centerZ) < RES / 10) {
				this->_density[_I(5, i, j)] = dens;
				this->_velX[_I(5, i, j)] = 2.0f;
			}
		}
	}

}

void Fluid::SimulateStep()
{
	GenerateSmoke();

	vel_step(N, vx_latest, vy_latest, vz_latest, vx_backup, vy_backup, vz_backup, visc, dt);
	dens_step(N, density_latest, density_backup, vx_latest, vy_latest, vz_latest, diff, dt);

}


void Fluid::Show()
{
	_renderer->FillTexture();
	_renderer->Render();
	
}


void Fluid::ClearBuffer(float* buf)
{
	for (int i = 0; i < SIZE; i++) {
		buf[i] = 0.0f;
	}
}

void Fluid::ClearSources(void)
{
	for (int i = 0; i < SIZE; i++) {
		sd[i] = su[i] = sv[i] = sw[i] = 0.0f;
	}
}

*/























