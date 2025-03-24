% Load the Java networking library
import java.net.*;

% Get the network interface of the computer
ni = NetworkInterface.getNetworkInterfaces();

% Initialize an empty list to store the MAC addresses
mac_list = {};

% Loop through the network interfaces
while ni.hasMoreElements()
    % Get the next network interface
    iface = ni.nextElement();
    
    % Get the hardware address of the interface
    mac = iface.getHardwareAddress();
    
    % If the hardware address is not empty, append it to the list
    if ~isempty(mac)
        mac_list{end+1} = dec2hex(mac)';
    end
end

% Print the list of MAC addresses
% disp(mac_list);

if ismember('36CFF66CF02F', mac_list) ||...
   ismember('F8E4E3B46C5C', mac_list) ||...
   ismember('FCAA14C8256D', mac_list) ||... % union liuyuan
   ismember('00FF964C5DED', mac_list) ||... % union jinziwei
   ismember('B4B5B6F798A7', mac_list) ||... % union guoyiwan
   ismember('34415DBF6B25', mac_list) ||... % union taoran
   ismember('088E907F0F98', mac_list) ||... % union qianlijuan
   ismember('00FF7FF46159', mac_list) ||... % union fanwenliang
   ismember('6C0B84AA21D0', mac_list) ||... % union fanwenliang
   ismember('00FFAB56CBB4', mac_list) ||... % union zhanglan
   ismember('00FF7916FB14', mac_list) ||... % union pengshuyi
   ismember('F4A475AE9ECF', mac_list) ||... % union zhangxinli B0-7B-25-5D-7D-9E
   ismember('7C8AE1A996E0', mac_list) ||... % union zhangxinli B0-7B-25-5D-7D-9E 
   ismember('1C1BB54BEE93', mac_list) ||... % union lujue
   ismember('0045E25D383B', mac_list) ||... % union zhanglan2
   ismember('6C0B840C1A99', mac_list) ||... % union liuyuan1
   ismember('1091D14EB5F6', mac_list) ||... % union liuyuan2
   ismember('04CF4BB1588C', mac_list) ||... % SP  laptop (win11)
   ismember('005056C00001', mac_list) ||... % Union wangjing
   ismember('E40D360B7A39', mac_list) ||... % Union wangjing   
   ismember('F8A2D6D3DD79', mac_list) ||... % Union wangjing
   ismember('6680726A1396', mac_list) ||... % ZHUyi
   ismember('6CF6DA379B55', mac_list)  % SP lenovo desktop (win11)
   disp('license checked');
else
   disp('The mac_list does not contain this PC');
   error('The mac_list must contain the MAC of this PC');
end