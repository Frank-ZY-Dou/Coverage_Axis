def read_VD(path):
    points = []
    radius = []
    with open(path, "r") as f:
        for line in f.readlines():
            line = line.strip('\n')  # 去掉列表中每一个元素的换行符
            line = line.split(' ')
            points.append([float(line[1]), float(line[2]), float(line[3])])
            radius.append([float(line[4])])
    return points, radius


def save_obj(path, verts, faces=None):
    verts = verts.tolist()

    with open(path, 'w') as f:
        for v in verts:
            f.write('v %f %f %f\n' %(v[0], v[1], v[2]))
        if faces is not None:
            faces = faces.tolist()
            for ff in faces:
                f.write('f %d %d %d\n' % (ff[0] + 1, ff[1] + 1, ff[2] + 1))
