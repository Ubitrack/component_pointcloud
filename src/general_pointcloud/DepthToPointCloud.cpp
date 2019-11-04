/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen, and individual
 * contributors as indicated by the @authors tag. See the
 * copyright.txt in the distribution for a full listing of individual
 * contributors.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 */

/**
 * @ingroup pcl_components
 * @file
 * Componenten to create a point cloud from a depth map
 *
 * 
 */

#include <utDataflow/TriggerComponent.h>
#include <utDataflow/TriggerInPort.h>
#include <utDataflow/TriggerOutPort.h>
#include <utDataflow/ComponentFactory.h>
#include <utMeasurement/Measurement.h>
#include <utVision/Image.h>
#include <utUtil/Logging.h>

static log4cpp::Category& logger(log4cpp::Category::getInstance("Ubitrack.PointCloud.DepthMaptoPointCloud"));

namespace Ubitrack { namespace Components {





        static bool transformation_project_internal(const Math::CameraIntrinsics<double> & calib,
                                                    const Math::Vector2d& xy,
                                                    Math::Vector2d& uv,
                                                    int *valid,
                                                    Math::Vector4d& J_xy)
        {

            float cx = -static_cast<float>(calib.matrix(0,2));
            float cy = -static_cast<float>(calib.matrix(1,2));
            float fx = static_cast<float>(calib.matrix(0,0));
            float fy = static_cast<float>(calib.matrix(1,1));

            float k1 = 0;
            float k2 =  0;
            float k3 =  0;
            float k4 =  0;
            float k5 =  0;
            float k6 =  0;
            switch (calib.radial_size) {
                case 6:
                    k6 = static_cast<float>(calib.radial_params[5]);
                case 5:
                    k5 = static_cast<float>(calib.radial_params[4]);
                case 4:
                    k4 = static_cast<float>(calib.radial_params[3]);
                case 3:
                    k3 = static_cast<float>(calib.radial_params[2]);
                case 2:
                    k2 = static_cast<float>(calib.radial_params[1]);
                    k1 = static_cast<float>(calib.radial_params[0]);
                default:
                    break;
            }

            float codx = 0.f; // center of distortion is set to 0 for Brown Conrady model
            float cody = 0.f;
            float p1 = static_cast<float>(calib.tangential_params[0]);
            float p2 = static_cast<float>(calib.tangential_params[1]);

            float max_radius_for_projection = 10.7f; // default value 1.7 corresponds to ~120 degree FoV, bigger to always make all with no check

            if (!(fx > 0.f && fy > 0.f))
            {
                LOG4CPP_ERROR(logger, "Expect both fx and fy are larger than 0, actual values are fx:" << (double)fx << ", fy: " << (double)fy);
                return false;
            }

            *valid = 1;

            float xp = xy[0] - codx;
            float yp = xy[1] - cody;

            float xp2 = xp * xp;
            float yp2 = yp * yp;
            float xyp = xp * yp;
            float rs = xp2 + yp2;
            if (rs > max_radius_for_projection * max_radius_for_projection)
            {
                *valid = 0;
                return true;
            }
            float rss = rs * rs;
            float rsc = rss * rs;
            float a = 1.f + k1 * rs + k2 * rss + k3 * rsc;
            float b = 1.f + k4 * rs + k5 * rss + k6 * rsc;
            float bi;
            if (b != 0.f)
            {
                bi = 1.f / b;
            }
            else
            {
                bi = 1.f;
            }
            float d = a * bi;

            float xp_d = xp * d;
            float yp_d = yp * d;

            float rs_2xp2 = rs + 2.f * xp2;
            float rs_2yp2 = rs + 2.f * yp2;

            xp_d += rs_2xp2 * p2 + 2.f * xyp * p1;
            yp_d += rs_2yp2 * p1 + 2.f * xyp * p2;

            float xp_d_cx = xp_d + codx;
            float yp_d_cy = yp_d + cody;

            uv[0] = xp_d_cx * fx + cx;
            uv[1] = yp_d_cy * fy + cy;

            // compute Jacobian matrix
            float dudrs = k1 + 2.f * k2 * rs + 3.f * k3 * rss;
            // compute d(b)/d(r^2)
            float dvdrs = k4 + 2.f * k5 * rs + 3.f * k6 * rss;
            float bis = bi * bi;
            float dddrs = (dudrs * b - a * dvdrs) * bis;

            float dddrs_2 = dddrs * 2.f;
            float xp_dddrs_2 = xp * dddrs_2;
            float yp_xp_dddrs_2 = yp * xp_dddrs_2;

            // compute d(u)/d(xp)
            J_xy[0] = fx * (d + xp * xp_dddrs_2 + 6.f * xp * p2 + 2.f * yp * p1);
            J_xy[1] = fx * (yp_xp_dddrs_2 + 2.f * yp * p2 + 2.f * xp * p1);
            J_xy[2] = fy * (yp_xp_dddrs_2 + 2.f * xp * p1 + 2.f * yp * p2);
            J_xy[3] = fy * (d + yp * yp * dddrs_2 + 6.f * yp * p1 + 2.f * xp * p2);

            return true;
        }

        static void invert_2x2(const Math::Vector4d& J, Math::Vector4d& Jinv)
        {
            float detJ = J[0] * J[3] - J[1] * J[2];
            float inv_detJ = 1.f / detJ;

            Jinv[0] = inv_detJ * J[3];
            Jinv[3] = inv_detJ * J[0];
            Jinv[1] = -inv_detJ * J[1];
            Jinv[2] = -inv_detJ * J[2];
        }

        static bool transformation_iterative_unproject(const Math::CameraIntrinsics<double> & calib,
                                                       const Math::Vector2d& uv,
                                                       Math::Vector2d& xy,
                                                       int *valid,
                                                       unsigned int max_passes)
        {
            *valid = 1;
            Math::Vector4d Jinv;
            Math::Vector2d best_xy{ 0.f, 0.f };
            float best_err = std::numeric_limits<float>::max();

            for (unsigned int pass = 0; pass < max_passes; pass++)
            {
                Math::Vector2d p;
                Math::Vector4d J;

                if (!(transformation_project_internal(calib, xy, p, valid, J)))
                {
                    return false;
                }
                if (*valid == 0)
                {
                    return true;
                }

                float err_x = uv[0] - p[0];
                float err_y = uv[1] - p[1];
                float err = err_x * err_x + err_y * err_y;
                if (err >= best_err)
                {
                    xy[0] = best_xy[0];
                    xy[1] = best_xy[1];
                    break;
                }

                best_err = err;
                best_xy[0] = xy[0];
                best_xy[1] = xy[1];
                invert_2x2(J, Jinv);
                if (pass + 1 == max_passes || best_err < 1e-22f)
                {
                    break;
                }

                float dx = Jinv[0] * err_x + Jinv[1] * err_y;
                float dy = Jinv[2] * err_x + Jinv[3] * err_y;

                xy[0] += dx;
                xy[1] += dy;
            }

            if (best_err > 1e-6f)
            {
                *valid = 0;
            }

            return true;
        }

        static bool transformation_unproject_internal(const Math::CameraIntrinsics<double> & calib,
                                                      const Math::Vector2d& uv,
                                                      Math::Vector2d& xy,
                                                      int *valid)
        {

            float cx = -static_cast<float>(calib.matrix(0,2));
            float cy = -static_cast<float>(calib.matrix(1,2));
            float fx = static_cast<float>(calib.matrix(0,0));
            float fy = static_cast<float>(calib.matrix(1,1));

            float k1 = 0;
            float k2 =  0;
            float k3 =  0;
            float k4 =  0;
            float k5 =  0;
            float k6 =  0;
            switch (calib.radial_size) {
                case 6:
                    k6 = static_cast<float>(calib.radial_params[5]);
                case 5:
                    k5 = static_cast<float>(calib.radial_params[4]);
                case 4:
                    k4 = static_cast<float>(calib.radial_params[3]);
                case 3:
                    k3 = static_cast<float>(calib.radial_params[2]);
                case 2:
                    k2 = static_cast<float>(calib.radial_params[1]);
                    k1 = static_cast<float>(calib.radial_params[0]);
                default:
                    break;
            }

            float codx = 0.f; // center of distortion is set to 0 for Brown Conrady model
            float cody = 0.f;
            float p1 = static_cast<float>(calib.tangential_params[0]);
            float p2 = static_cast<float>(calib.tangential_params[1]);

            if (!(fx > 0.f && fy > 0.f))
            {
                LOG4CPP_ERROR(logger, "Expect both fx and fy are larger than 0, actual values are fx:" << (double)fx << ", fy: " << (double)fy);
                return false;
            }

            // correction for radial distortion
            float xp_d = (uv[0] - cx) / fx - codx;
            float yp_d = (uv[1] - cy) / fy - cody;

            float rs = xp_d * xp_d + yp_d * yp_d;
            float rss = rs * rs;
            float rsc = rss * rs;
            float a = 1.f + k1 * rs + k2 * rss + k3 * rsc;
            float b = 1.f + k4 * rs + k5 * rss + k6 * rsc;
            float ai;
            if (a != 0.f)
            {
                ai = 1.f / a;
            }
            else
            {
                ai = 1.f;
            }
            float di = ai * b;

            xy[0] = xp_d * di;
            xy[1] = yp_d * di;

            // approximate correction for tangential params
            float two_xy = 2.f * xy[0] * xy[1];
            float xx = xy[0] * xy[0];
            float yy = xy[1] * xy[1];

            xy[0] -= (yy + 3.f * xx) * p2 + two_xy * p1;
            xy[1] -= (xx + 3.f * xx) * p1 + two_xy * p2;

            // add on center of distortion
            xy[0] += codx;
            xy[1] += cody;

            return transformation_iterative_unproject(calib, uv, xy, valid, 20);
        }

        bool transformation_unproject(const Math::CameraIntrinsics<double> & calib,
                                      const Math::Vector2d& point2d,
                                      const float depth,
                                      Math::Vector3d& point3d,
                                      int *valid)
        {
            if (depth == 0.f)
            {
                point3d[0] = 0.f;
                point3d[1] = 0.f;
                point3d[2] = 0.f;
                *valid = 0;
                return true;
            }

            Math::Vector2d xy{ point3d[0], point3d[1] };
            if (!transformation_unproject_internal(calib, point2d, xy, valid))
            {
                return false;
            }

            point3d[0] = xy[0] * depth;
            point3d[1] = xy[1] * depth;
            point3d[2] = depth;

            return true;
        }

        bool transformation_project(const Math::CameraIntrinsics<double> & calib,
                                    const Math::Vector3d& point3d,
                                    Math::Vector2d& point2d,
                                    int *valid)
        {
            if (point3d[2] <= 0.f)
            {
                point2d[0] = 0.f;
                point2d[1] = 0.f;
                *valid = 0;
                return true;
            }

            Math::Vector2d xy;
            xy[0] = point3d[0] / point3d[2];
            xy[1] = point3d[1] / point3d[2];

            Math::Vector4d _J; // unused
            return transformation_project_internal(calib, xy, point2d, valid, _J);
        }


/**
 * @ingroup pcl_components
 */
class PCL_DepthToPointCloud
	: public Dataflow::TriggerComponent
{
public:
	/**
	 * UTQL component constructor.
	 *
	 * @param sName Unique name of the component.
	 * @param subgraph UTQL subgraph
	 */
    PCL_DepthToPointCloud( const std::string& sName, boost::shared_ptr< Graph::UTQLSubgraph > pConfig )
		: Dataflow::TriggerComponent( sName, pConfig )
		, m_inPortImage( "DepthImage", *this )
		, m_inPortIntrinsics( "Intrinsics", *this )
		, m_outPort( "Output", *this )
    {
    }

	/** Method that computes the result. */
	void compute( Measurement::Timestamp t )
	{
	    if(m_xy_table.size() == 0) {
            create_xy_table(*m_inPortIntrinsics.get());
            LOG4CPP_INFO(logger, "initializing xy table, example vectors");
            LOG4CPP_INFO(logger,m_xy_table[0]);
	    }

	    cv::Mat depthImage = m_inPortImage.get()->Mat();

	    std::vector<Math::Vector3d> result;

	    uint16_t* depthData = reinterpret_cast<uint16_t*>( depthImage.data );

	    double depthScale = 0.001;

	    for(size_t i=0;i<m_pointCount;i++) {
	        uint16_t depthInt = depthData[i];
	        double depth = static_cast<double>(depthInt);
	        Math::Vector2d pointVec = m_xy_table[i];

	        if(depth == 0)
                continue;
            double realDepth = depthScale * depth;
            Math::Vector3d xyz = Math::Vector3d(pointVec(0) * realDepth, pointVec(1) * realDepth, -realDepth);
            result.push_back(xyz);


	    }





		m_outPort.send( Measurement::PositionList( t, result ) );
	}

protected:
	/** Input port A of the component. */
	Dataflow::TriggerInPort< Measurement::ImageMeasurement > m_inPortImage;

	/** Input port B of the component. */
	Dataflow::TriggerInPort< Measurement::CameraIntrinsics > m_inPortIntrinsics;

	/** Output port of the component. */
	Dataflow::TriggerOutPort< Measurement::PositionList > m_outPort;

	std::vector<Math::Vector2d >  m_xy_table;
	size_t m_pointCount=0;


    void create_xy_table(const Math::CameraIntrinsics<double> & calib)
    {
        int width = calib.dimension[0];
        int height = calib.dimension[1];
        m_pointCount= width*height;

        m_xy_table.clear();
        m_xy_table.reserve(m_pointCount);



        Math::Vector2d p;
        Math::Vector3d ray;
        int valid;

        // precompute xy lookup table (Vec2(x,y,1) are vectors pointing to the image plane at distance 1 unit)
        for (int y = 0, idx = 0; y < height; y++)
        {
            p(1) = (float)y;
            for (int x = 0; x < width; x++, idx++)
            {
                p(0) = (float)x;

                if (transformation_unproject(calib, p, 1.f, ray, &valid)) {
                    if (valid)
                    {
                        m_xy_table.push_back(Math::Vector2d(ray(0), -ray(1)));
                    }
                    else
                    {
                        m_xy_table.push_back(Math::Vector2d(0,0));
                    }
                }
            }
        }

    }
};


UBITRACK_REGISTER_COMPONENT( Dataflow::ComponentFactory* const cf ) {
	cf->registerComponent< PCL_DepthToPointCloud > ( "PCL_DepthToPointCloud" );


}

} } // namespace Ubitrack::Components
