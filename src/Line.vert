#version 410 core
in vec3 Position;
in mat4 InstanceMatrix;

uniform mat4 Projection;
uniform mat4 View;
uniform vec4 Color;

out vec4 fColor;

void main(void)
{
	fColor = Color;
	gl_Position = Projection * View * InstanceMatrix * vec4(Position, 1.0);
}
