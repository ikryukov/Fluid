#version 410 core
in vec3 Position;
in vec4 InstanceColor;
in mat4 InstanceMatrix;

uniform mat4 Projection;
uniform mat4 View;

out vec4 fColor;

void main(void)
{
	fColor = InstanceColor;
	gl_Position = Projection * View * InstanceMatrix * vec4(Position, 1.0);
}
