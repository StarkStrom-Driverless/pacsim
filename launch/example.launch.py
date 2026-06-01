import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.substitutions import Command, LaunchConfiguration
from launch.event_handlers import OnProcessExit
from launch.actions import RegisterEventHandler, LogInfo, EmitEvent, DeclareLaunchArgument
from launch.events import Shutdown

def getFullFilePath(name, dir):
    ret = os.path.join(
        get_package_share_directory('pacsim'),
        dir,
        name
    )
    return ret


def generate_launch_description():
  track_name_arg = DeclareLaunchArgument(
      'track_name',
      default_value=getFullFilePath("FSG23.yaml", "tracks"),
      description='Full path to the track YAML file'
  )
  discipline_arg = DeclareLaunchArgument(
      'discipline',
      default_value='autocross',
      description='Driving discipline'
  )
  headless_arg = DeclareLaunchArgument(
      'headless',
      default_value='false',
      description='Run in headless mode (no RViz)'
  )
  report_file_dir_arg = DeclareLaunchArgument(
      'report_file_dir',
      default_value='/tmp',
      description='Directory to write the evaluation report file'
  )

  track_name = LaunchConfiguration('track_name')
  discipline = LaunchConfiguration('discipline')
  report_file_dir = LaunchConfiguration('report_file_dir')

  track_frame = "map"
  realtime_ratio = 1.0

  xacro_file_name = 'separate_model.xacro'
  xacro_path = getFullFilePath(xacro_file_name, "urdf")

  nodePacsim = Node(
          package='pacsim',
          namespace='pacsim',
          executable='pacsim_node',
          name='pacsim_node',
          parameters = [{'use_sim_time':True}, {"track_name" : track_name}, {"grip_map_path" : getFullFilePath("gripMap.yaml", "tracks")}, {"track_frame" : track_frame}, {"realtime_ratio" : realtime_ratio},
                        {"report_file_dir" : report_file_dir}, {"main_config_path" : getFullFilePath("mainConfig.yaml", dir="config")}, {"perception_config_path" : getFullFilePath("perception.yaml", dir="config")},
                        {"sensors_config_path" : getFullFilePath("sensors.yaml", dir="config")}, {"vehicle_model_config_path" : getFullFilePath("vehicleModel.yaml", dir="config")}, {"discipline" : discipline}],
          output="screen",
          emulate_tty=True)
  nodePacsimShutdownEventHandler = RegisterEventHandler(
        OnProcessExit(
            target_action=nodePacsim,
            on_exit=[
                LogInfo(msg=('Pacsim closed')),
                EmitEvent(event=Shutdown(
                    reason='Pacsim closed, shutdown whole launchfile')),
            ]
        )
  )
  robot_state_publisher= Node(
          package='robot_state_publisher',
          executable='robot_state_publisher',
          name='robot_state_publisher',
          output='screen',
          parameters=[{'use_sim_time':True}, {'publish_frequency':float(1000),
                      'robot_description': Command(['xacro',' ',xacro_path])
                      }],
          arguments=[xacro_path])
  return LaunchDescription([track_name_arg, discipline_arg, headless_arg, report_file_dir_arg, nodePacsim, nodePacsimShutdownEventHandler, robot_state_publisher])
